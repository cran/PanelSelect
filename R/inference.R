# 1. Coordinate scales -------------------------------------------------------

# Root mean square of the finite nonzero entries of each design column, named
# by coefficient: selection columns over all rows, outcome columns over the
# selected rows. It carries the units of each coefficient's regressor.
design_scales = function(prepared) {
  rms = function(matrix) vapply(seq_len(ncol(matrix)), function(j) {
    value = matrix[, j]; value = value[is.finite(value) & value != 0]
    if (length(value)) sqrt(mean(value^2)) else 1
  }, numeric(1))
  c(if (prepared$has_selection && ncol(prepared$w))
      stats::setNames(rms(prepared$w), paste0("selection.", colnames(prepared$w))),
    if (ncol(prepared$x)) stats::setNames(rms(prepared$x[prepared$selected, , drop = FALSE]),
      paste0("outcome.", colnames(prepared$x))))
}

# Floors of the split and staged integration score budgets, tau*(floor +
# |score|), named by parameter: the design scale for regression coefficients
# and one for log-SD and atanh-correlation coordinates. Rescaling a regressor
# rescales its scores, errors and floor together, so refinement decisions do
# not depend on covariate units. Floors are cached in an evaluator's private
# workspace, whose data never change.
score_floors = function(prepared, theta, workspace = NULL) {
  if (is.environment(workspace) && !is.null(workspace$score_floors)) return(workspace$score_floors)
  floors = stats::setNames(rep(1, length(theta)), names(theta))
  scales = if (!is.null(prepared$x)) design_scales(prepared)
  coefficient = intersect(names(theta), names(scales))
  floors[coefficient] = scales[coefficient]
  if (is.environment(workspace)) workspace$score_floors = floors
  floors
}

# 2. Observed information ----------------------------------------------------

# Score-difference steps at one endpoint, in free optimizer coordinates.
# Gaussian outcome coefficients move their linear predictor by about 1e-3
# residual SDs; other coefficients move their index by about 1e-3. Dividing
# by the design scale respects regressor units. Log-SD and atanh-correlation
# steps are unit free. Keep these endpoint scales fixed throughout the stencil.
information_steps = function(info, prepared, theta, scale = design_scales(prepared)) {
  free = info$name[!info$fixed]
  step = stats::setNames(rep(1e-3, length(free)), free)
  coefficient = intersect(free, names(scale))
  step[coefficient] = 1e-3/scale[coefficient]
  if (prepared$family == "gaussian") {
    outcome = coefficient[startsWith(coefficient, "outcome.")]
    step[outcome] = step[outcome]*theta[["sd_error"]]
  }
  step
}

information_matrix = function(eta, evaluator, on_error = NULL) {
  if (!length(eta)) return(matrix(numeric(), 0L, 0L))
  objective = function(value) -sum(evaluator$evaluate(value, trial = FALSE)$loglik)
  gated = !is.null(evaluator$information_gate_tolerance)
  stencil = new.env(parent = emptyenv())
  stencil_points = new.env(parent = emptyenv())
  coordinate_error = 0
  coordinate_matches = list()
  key = function(value) paste(as.integer(serialize(as.numeric(value), NULL, version = 2)), collapse = ',')
  gradient = function(value) {
    result = evaluator$evaluate(value, trial = FALSE)
    if (gated) {
      stencil[[key(value)]] = -colSums(result$reference_optimizer_score)
      stencil_points[[key(value)]] = as.numeric(value)
    }
    -colSums(result$optimizer_score)
  }
  tryCatch({
    steps = unname(evaluator$information_steps(eta))
    if (!is.numeric(steps) || is.complex(steps) || !is.null(dim(steps)) ||
        length(steps) != length(eta))
      stop("Information steps do not match the free parameters", call. = FALSE)
    if (any(!is.finite(steps) | steps <= 0))
      stop("Information steps must be positive and finite", call. = FALSE)
    # A rounded-away or overflowing perturbation cannot identify curvature.
    # Reject it through the ordinary unavailable-information path, without a
    # unit-dependent floor or a silently enlarged step. Check the GL half-step
    # as well as the native stencil before evaluating either one.
    for (step in if (gated) list(steps, steps/2) else list(steps)) {
      plus = eta + step; minus = eta - step
      if (any(!is.finite(c(eta, plus, minus))) || any(plus == eta | minus == eta))
        stop("Observed-information perturbations are not representable at this endpoint", call. = FALSE)
    }
    if (!is.null(evaluator$anchor)) {
      restore = evaluator$anchor(eta)
      on.exit(restore(), add = TRUE)
    }
    information = stats::optimHess(eta, objective, gradient, control = list(ndeps = steps))
    if (!is.null(evaluator$information_gate_tolerance)) {
      refined = stats::optimHess(eta, objective, gradient, control = list(ndeps = steps/2))
      lower_information = function(step) {
        answer = matrix(0, length(eta), length(eta), dimnames = list(names(eta), names(eta)))
        for (j in seq_along(eta)) {
          plus = minus = eta; plus[j] = plus[j] + step[j]; minus[j] = minus[j] - step[j]
          reference = function(point) {
            exact = stencil[[key(point)]]
            if (!is.null(exact)) return(exact)
            # optimHess restores coordinates by sequential floating-point
            # additions/subtractions. A nonvaried coordinate can differ from
            # eta by a few ulps. Match only such roundoff, never a nearby step.
            keys = ls(stencil_points, all.names = TRUE)
            distance = vapply(keys, function(name) {
              actual = stencil_points[[name]]
              max(abs(actual - point)/pmax(1, abs(actual), abs(point)))
            }, numeric(1))
            nearest = which(is.finite(distance) & distance <= 16*.Machine$double.eps)
            if (length(nearest) != 1L)
              integration_error("Observed-information accuracy check: reference stencil was not retained within roundoff")
            coordinate_error <<- max(coordinate_error, distance[nearest])
            coordinate_matches[[length(coordinate_matches) + 1L]] <<- list(requested = as.numeric(point),
              actual = stencil_points[[keys[nearest]]],
              difference = stencil_points[[keys[nearest]]] - as.numeric(point),
              scaled_error = distance[nearest])
            stencil[[keys[nearest]]]
          }
          answer[, j] = (reference(plus) - reference(minus))/(2*step[j])
        }
        (answer + t(answer))/2
      }
      accuracy = tryCatch(information_accuracy_gate(information, refined,
        lower_information(steps), lower_information(steps/2),
        evaluator$information_gate_tolerance, steps),
        panelselect_integration_error = function(error) {
          if (isTRUE(evaluator$staged)) stop(staged_refinement_condition(conditionMessage(error),
            diagnostics = list(class = class(error), cause = error[['diagnostics']], phase = 'information')))
          stop(error)
        })
      if (!is.null(evaluator$information_gate_basis)) accuracy$basis = evaluator$information_gate_basis
      accuracy$stencil_coordinate_max_scaled_roundoff = coordinate_error
      accuracy$stencil_coordinate_tolerance = 16*.Machine$double.eps
      accuracy$stencil_coordinate_matches = coordinate_matches
      # The attribute keeps its serialized name for saved fits.
      attr(information, "checked_accuracy") = accuracy
    }
    information
  }, error = function(error) {
    # One exiting handler preserves the typed request. Re-throwing from an
    # earlier handler in this tryCatch would be caught by its generic handler.
    if (inherits(error, 'panelselect_mesh_required')) stop(error)
    if (!is.null(on_error)) on_error(conditionMessage(error))
    warning("Observed information is unavailable: ", conditionMessage(error), call. = FALSE)
    matrix(NA_real_, length(eta), length(eta), dimnames = list(names(eta), names(eta)))
  })
}

# Step/order and inverse-stability gate for split and staged observed information.
information_accuracy_gate = function(information, refined, lower, lower_refined, tolerance,
                                    steps = rep(.001, nrow(information))) {
  if (!identical(dim(information), dim(refined)) || !identical(dim(information), dim(lower)) ||
      !identical(dim(information), dim(lower_refined)) ||
      any(!is.finite(c(information, refined, lower, lower_refined))))
    integration_error('Observed-information accuracy check: information matrices are unavailable or not finite')
  # Entry changes standardized by the diagonals, |a_jk - b_jk|/(d_j d_k) with
  # d_j = sqrt(max(|a_jj|, |b_jj|)), so rescaled optimizer coordinates leave
  # the statistics unchanged. A coordinate without curvature in either matrix
  # (d_j = 0) is compared absolutely.
  scaled_change = function(a, b) {
    d = sqrt(pmax(abs(diag(a)), abs(diag(b)))); d[d == 0] = 1
    max(c(0, abs(a - b)/outer(d, d)))
  }
  discrepancy = scaled_change(information, refined)
  order_change = max(scaled_change(information, lower), scaled_change(refined, lower_refined))
  if (max(discrepancy, order_change) > tolerance)
    integration_error(paste('Observed-information accuracy check: step/order refinement unresolved:',
      signif(discrepancy, 5), '/', signif(order_change, 5)))
  inverse = function(x) tryCatch({
    value = chol2inv(chol(x))
    if (any(!is.finite(value)) || any(diag(value) <= 0)) return(NULL)
    value
  }, error = function(e) NULL)
  native_covariance = inverse(information)
  covariance_check = list(status = 'UNAVAILABLE_NATIVE_NOT_PD',
    qualification = 'An indefinite/singular native information matrix remains available for derivative diagnostics; observed covariance remains unavailable.')
  status = 'DERIVATIVE_REFINEMENT_PASS_COVARIANCE_UNAVAILABLE'
  if (!is.null(native_covariance)) {
    refinements = list(half_step = refined, predecessor = lower,
      half_step_predecessor = lower_refined)
    native_sd = sqrt(diag(native_covariance)); scale = outer(native_sd, native_sd)
    if (any(!is.finite(scale)) || any(scale <= 0))
      integration_error('Observed-information accuracy check: covariance scaling is not finite')
    comparisons = lapply(names(refinements), function(name) {
      alternative = inverse(refinements[[name]])
      if (is.null(alternative))
        integration_error(paste('Observed-information accuracy check: refined information is not positive definite:', name))
      se_change = max(abs(sqrt(diag(alternative))/native_sd - 1))
      covariance_change = max(abs(alternative - native_covariance)/scale)
      if (!is.finite(se_change) || !is.finite(covariance_change) ||
          max(se_change, covariance_change) > tolerance)
        integration_error(paste('Observed-information accuracy check: covariance refinement unresolved:', name,
          'relative SE', signif(se_change, 5), 'standardized covariance', signif(covariance_change, 5)))
      list(rule = name, maximum_relative_SE_change = se_change,
        maximum_standardized_covariance_change = covariance_change)
    })
    covariance_check = list(status = 'PASS', tolerance = tolerance, comparisons = comparisons,
      qualification = 'Observed-information inverse/SE stability at this endpoint; explicit OPG covariance has a separate score-inverse stability question.')
    status = 'PASS'
  }
  list(status = status, scaled_max_change = discrepancy, scaled_order_change = order_change, tolerance = tolerance,
    observed_covariance_refinement = covariance_check,
    native_optimizer_step = steps, refined_optimizer_step = steps/2,
    canonical_information = information, predecessor_information = lower,
    half_step_information = refined, half_step_predecessor_information = lower_refined,
    basis = 'Complete-value numerical scores on both information stencils at canonical and predecessor quadrature orders; covariance remains based on the original canonical native stencil.')
}

# 3. Covariance estimators and coordinate transformation ----------------------

positive_inverse = function(matrix) {
  if (!nrow(matrix)) return(matrix)
  if (any(!is.finite(matrix))) return(matrix * NA_real_)
  tryCatch({
    value = chol2inv(chol(matrix))
    dimnames(value) = dimnames(matrix)
    value
  }, error = function(error) matrix * NA_real_)
}

covariance_matrices = function(theta, info, score, information) {
  free = !info$fixed
  jacobian = parameter_jacobian(theta, info)[free]
  optimizer_score = sweep(score[, free, drop = FALSE], 2L, jacobian, `*`)
  meat = crossprod(optimizer_score)
  observed = positive_inverse(information)
  optimizer = list(observed = observed, opg = positive_inverse(meat),
                    sandwich = observed %*% meat %*% observed)
  natural = lapply(optimizer, function(value) {
    result = matrix(0, nrow(info), nrow(info), dimnames = list(info$name, info$name))
    result[free, free] = value * outer(jacobian, jacobian)
    result
  })
  list(natural = natural, optimizer = optimizer)
}

# Observed information at the fitted optimizer coordinates, computed exactly as
# a vcov = "observed" fit computes it. A staged fit whose retained mesh cannot
# serve the information stencil warns and returns an unavailable (NA) matrix;
# the reason is returned as `error`, as a fit records information_error.
observed_information = function(object) {
  optimizer = object[["optimizer"]]
  error = NULL
  evaluator = fit_evaluator(object$coefficients, object$parameter_info, object$prepared,
                             object$control, object$fixed)
  information = tryCatch({
    if (!is.null(object$control$integration_mesh))
      validate_staged_integration(object$coefficients, object$prepared, object$control,
        object$fixed, mesh = object$control$integration_mesh, information = TRUE)
    information_matrix(optimizer[["par"]], evaluator, on_error = function(message) error <<- message)
  }, panelselect_mesh_required = function(condition) {
    error <<- conditionMessage(condition)
    warning('Observed information requires a refined integration mesh; use update/refit with vcov="observed": ',
      conditionMessage(condition), call. = FALSE)
    matrix(NA_real_, length(optimizer[['par']]), length(optimizer[['par']]),
      dimnames = list(names(optimizer[['par']]), names(optimizer[['par']])))
  })
  list(information = information, error = error)
}

covariance_pair = function(object, type, compute = TRUE) {
  type = match.arg(type, c("observed", "opg", "sandwich"))
  optimizer = object[["optimizer"]]
  pair = list(natural = object[["covariance"]][[type]],
    optimizer = optimizer[["covariance"]][[type]], eta = optimizer[["par"]])
  if (compute && type != "opg" && identical(optimizer[["information_computed"]], FALSE)) {
    information = observed_information(object)$information
    covariance = covariance_matrices(object$coefficients, object$parameter_info,
                                     object$score, information)
    pair$natural = covariance$natural[[type]]
    pair$optimizer = covariance$optimizer[[type]]
  }
  pair
}

# 4. Store and reuse observed information ------------------------------------

#' Compute and store observed information for a fitted model
#'
#' Augment an OPG fit with the observed Hessian and its covariance estimators,
#' without refitting. Assign the returned object to retain the calculation for
#' later inference and stationarity diagnostics.
#'
#' @param object A fitted model of class `panelselect`, returned by
#'   [panel_select()]. Its fitted parameters, scores, controls and data define
#'   the information calculation.
#'
#' @details
#' The default `vcov = "opg"` fit does not compute a post-fit observed Hessian.
#' `add_information()` differentiates the full score at the fitted free
#' optimizer coordinates to obtain the negative log-likelihood Hessian. This is
#' the observed information of the fitted numerical objective; it does not
#' certify quadrature accuracy, a global maximum or calibrated intervals.
#'
#' The returned fit stores the information, any integration-accuracy record,
#' its computed/unavailable status and error reason, and the observed, OPG and
#' sandwich covariances on both optimizer and natural parameter scales. Its
#' coefficients, likelihood and default covariance choice are unchanged.
#' [check_stationarity()] can then report `Q_Hessian` and Newton diagnostics
#' when the stored information and its inverse are finite and positive
#' definite. Missing or invalid curvature remains unavailable; no OPG or
#' sandwich matrix substitutes for the observed inverse.
#'
#' Explicit observed/sandwich requests in [vcov()], [confint()], [summary()] or
#' predictions can calculate missing information for that request. They do not
#' store it in the caller's original object. Use `fit = add_information(fit)`
#' once to reuse the stored calculation in subsequent requests.
#'
#' A fit whose information is already marked as computed is returned unchanged,
#' including when its stored information is unavailable. This also applies to
#' information stored by older versions, including before the Gaussian
#' outcome-unit correction. Use `update(fit, vcov = "observed")` to refit and
#' recompute information with current numerical steps.
#'
#' A staged GL fit whose retained mesh cannot serve the information stencil
#' warns and stores an unavailable result. It does not silently refine the mesh
#' or replace the estimator. An explicit refit with `vcov = "observed"` prepares
#' a mesh for the information calculation, subject to the requested accuracy
#' and resource limits.
#'
#' @return A `panelselect` fit with stored information and covariance results,
#'   or the unchanged fit when information was already computed. Assign this
#'   return value to preserve the added information.
#' @family diagnostic functions
#' @seealso [panelselect-methods]
#' @examples
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#' fit = panel_select(y ~ x, observed ~ x + z, dat, "person", "linear")
#' fit = add_information(fit)
#' sqrt(diag(vcov(fit, "observed")))
#' check_stationarity(fit)$Q_Hessian
#' @export
add_information = function(object) {
  if (!inherits(object, "panelselect")) stop("object must be a panelselect fit", call. = FALSE)
  object = current_fit(object)
  if (!identical(object$optimizer[["information_computed"]], FALSE)) return(object)
  computed = observed_information(object)
  covariance = covariance_matrices(object$coefficients, object$parameter_info, object$score,
                                   computed$information)
  object$optimizer$information = computed$information
  object$optimizer$information_computed = TRUE
  object$optimizer["information_error"] = list(computed$error)
  object$optimizer$covariance = covariance$optimizer
  object$covariance = covariance$natural
  diagnostic = check_stationarity(object)
  object$optimizer$stationarity = diagnostic$opg
  object$optimizer$stationarity$Q_Hessian = diagnostic$Q_Hessian
  object
}

# 5. Coefficient uncertainty on natural and optimizer scales -----------------

# A covariance request is an estimator name or a named natural-scale matrix,
# validated as for predictions. A supplied matrix has no optimizer coordinates,
# so its correlation and SD intervals use the delta-method natural SEs.
covariance_request = function(object, vcov) {
  if (is.matrix(vcov)) return(list(natural = prediction_covariance(object, vcov),
    optimizer = NULL, eta = NULL))
  covariance_pair(object, vcov)
}

coefficient_scales = function(object, covariance) {
  estimate = object$coefficients; info = object$parameter_info
  V = covariance$natural; parameter_names = names(estimate); p = length(estimate)
  named = !anyDuplicated(parameter_names) && !anyNA(parameter_names) &&
    length(parameter_names) == p && !anyDuplicated(rownames(V)) && !anyDuplicated(colnames(V)) &&
    setequal(rownames(V), parameter_names) && setequal(colnames(V), parameter_names)
  if (!is.numeric(estimate) || is.complex(estimate) || !is.matrix(V) || !is.numeric(V) ||
      is.complex(V) || !identical(dim(V), c(p,p)) || !named)
    stop("Natural covariance must be real and have unique row/column names matching coefficients", call. = FALSE)
  se = sqrt(diag(V[parameter_names, parameter_names, drop = FALSE])); se[info$fixed] = NA_real_
  jacobian = parameter_jacobian(estimate, info)
  eta = estimate
  eta[info$type == "sd"] = log(estimate[info$type == "sd"])
  eta[info$type == "correlation"] = atanh(estimate[info$type == "correlation"])
  eta_se = se/jacobian
  free = !info$fixed; free_names = names(estimate)[free]
  retained = covariance$eta; V = covariance$optimizer
  aligned = length(free_names) > 0L && is.numeric(retained) && !is.complex(retained) &&
    length(retained) == length(free_names) &&
    !anyDuplicated(names(retained)) && setequal(names(retained), free_names) &&
    is.matrix(V) && is.numeric(V) && !is.complex(V) && identical(dim(V), c(length(free_names), length(free_names))) &&
    !anyDuplicated(rownames(V)) && !anyDuplicated(colnames(V)) &&
    setequal(rownames(V), free_names) && setequal(colnames(V), free_names)
  if (aligned) {
    retained = retained[free_names]
    forward = parameter_from_optimizer(retained, estimate, info)[free]
    aligned = all(is.finite(retained)) && identical(unname(forward), unname(estimate[free]))
  }
  if (aligned) {
    eta[free] = retained
    eta_se[free] = sqrt(diag(V[free_names, free_names, drop = FALSE]))
    se[free] = abs(jacobian[free])*eta_se[free]
  }
  list(estimate = estimate, se = se, eta = eta, eta_se = eta_se)
}

# `covariance` is an estimator name, a natural-scale matrix, or a computed
# covariance_request() result, so one method call computes it only once.
coefficient_table = function(object, covariance, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie strictly between zero and one", call. = FALSE)
  if (!is.list(covariance)) covariance = covariance_request(object, covariance)
  scales = coefficient_scales(object, covariance)
  estimate = scales$estimate; se = scales$se
  info = object$parameter_info
  critical = stats::qnorm((1 + level) / 2)
  z = estimate / se
  lower = estimate - critical * se
  upper = estimate + critical * se
  sd = info$type == "sd" & !info$fixed
  correlation = info$type == "correlation" & !info$fixed
  lower[sd] = exp(scales$eta[sd] - critical * scales$eta_se[sd])
  upper[sd] = exp(scales$eta[sd] + critical * scales$eta_se[sd])
  z[sd] = NA_real_
  lower[correlation] = tanh(scales$eta[correlation] - critical * scales$eta_se[correlation])
  upper[correlation] = tanh(scales$eta[correlation] + critical * scales$eta_se[correlation])
  z[correlation] = scales$eta[correlation] / scales$eta_se[correlation]
  cbind(estimate = estimate, se = se, z = z, p = 2 * stats::pnorm(-abs(z)),
        lci = lower, uci = upper)
}
