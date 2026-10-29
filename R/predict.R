# 1. Model targets and design inputs

# Shared prediction inputs; actual response/joint equations belong to each model.
model_predictions = function(theta, x, w, offset_x, offset_w, family,
                              has_selection, u = NULL, v = NULL, du = NULL, dv = NULL,
                              targets = c('response', 'selection', 'joint', 'selected',
                                          'link', 'selection_link'),
                              structure = 'panel', model_id = NULL) {
  model = resolve_model(family, has_selection, structure, model_id)
  if (structure == 'cross_sectional' && (!is.null(u) || !is.null(v)))
    stop('Cross-sectional models do not support history or random effects', call. = FALSE)
  n = nrow(x); p = length(theta)
  zero = matrix(0, n, p, dimnames = list(NULL, names(theta)))
  derivative = function(name) {
    g = zero
    if (name %in% names(theta)) g[, name] = 1
    g
  }
  value = function(name) if (name %in% names(theta)) unname(theta[name]) else 0
  gradient_b = zero; gradient_a = zero
  bx = if (ncol(x)) paste0("outcome.", colnames(x)) else character()
  bw = if (ncol(w)) paste0("selection.", colnames(w)) else character()
  if (ncol(x)) gradient_b[, bx] = x
  if (ncol(w)) gradient_a[, bw] = w
  b = (if (ncol(x)) drop(x %*% theta[bx]) else rep(0, n)) + offset_x
  a = (if (has_selection && ncol(w)) drop(w %*% theta[bw]) else rep(0, n)) + offset_w
  delta = value("sd_re_selection"); lambda = value("sd_re_outcome")
  error = value("sd_error"); rho = value("cor_re"); tau = value("cor_error")
  gradient_delta = derivative("sd_re_selection"); gradient_lambda = derivative("sd_re_outcome")
  gradient_error = derivative("sd_error"); gradient_rho = derivative("cor_re")
  gradient_tau = derivative("cor_error")
  link = list(fit = b, gradient = gradient_b)
  selection_link = list(fit = a, gradient = gradient_a)
  if (all(targets %in% c('link', 'selection_link')))
    return(list(link = link, selection_link = selection_link)[targets])
  conditional = !is.null(v)
  # At explicit history nodes a/b become A/B. Population selection instead uses
  # its marginal standardized index; the gradient matrices follow the same map.
  if (conditional) {
    if (is.null(du)) du = zero
    if (is.null(dv)) dv = zero
    a = a + delta*u; b = b + lambda*v
    gradient_a = gradient_a + gradient_delta*u + delta*du
    gradient_b = gradient_b + gradient_lambda*v + lambda*dv
    selection_scale = 1
  } else {
    selection_scale = sqrt(1 + delta^2)
    a = a/selection_scale
    gradient_a = gradient_a/selection_scale - gradient_delta*(a*delta/selection_scale^2)
  }
  logP = if (has_selection) stats::pnorm(a, log.p = TRUE) else rep(0, n)
  P = exp(logP)
  logGP = if (has_selection) gradient_a*prediction_mills(a) else zero
  GP = logGP*P
  selection = list(fit = P, gradient = GP, log_fit = logP,
                    log_gradient = logGP)
  if (!any(targets %in% c('response', 'joint', 'selected')))
    return(list(selection = selection, link = link, selection_link = selection_link)[targets])
  joint_needed = any(targets %in% c('joint', 'selected'))
  context = list(a = a, b = b, gradient_b = gradient_b, gradient_a = gradient_a, zero = zero,
    delta = delta, lambda = lambda, error = error, rho = rho, tau = tau,
    gradient_delta = gradient_delta, gradient_lambda = gradient_lambda,
    gradient_error = gradient_error, gradient_rho = gradient_rho,
    gradient_tau = gradient_tau, conditional = conditional,
    selection_scale = selection_scale, selection = selection, joint_needed = joint_needed)
  result = model$predictions(context)
  if (!has_selection) result$joint = result$selected = result$response
  list(response = result$response, selection = selection, joint = result$joint,
       selected = result$selected, link = link, selection_link = selection_link)[targets]
}

# Prediction uses natural parameters and preserves the fitted design schema.
prediction_design = function(design, data) {
  if (is.null(design)) return(list(matrix = matrix(0, nrow(data), 0L),
                                 offset = rep(0, nrow(data))))
  terms = stats::delete.response(design$terms)
  frame = stats::model.frame(terms, strip_contrasts(data), na.action = stats::na.pass,
                             xlev = design$xlevels, drop.unused.levels = FALSE)
  x = stats::model.matrix(terms, frame, contrasts.arg = design$contrasts)
  if (!identical(colnames(x), design$columns))
    stop("Prediction design does not match the fitted columns", call. = FALSE)
  offset = stats::model.offset(frame)
  if (any(is.infinite(x)) || any(is.infinite(offset)))
    stop("Prediction covariates and offsets must not be infinite", call. = FALSE)
  list(matrix = x, offset = if (is.null(offset)) rep(0, nrow(x)) else offset)
}

prediction_mills = function(x) {
  answer = rep(NA_real_, length(x))
  good = is.finite(x)
  answer[good] = normal_mills(x[good])
  answer
}

prediction_bvn = function(a, b, r, da, db, dr) {
  n = length(a); p = ncol(da)
  logp = rep(NA_real_, n); score = matrix(NA_real_, n, p, dimnames = dimnames(da))
  good = which(is.finite(a) & is.finite(b) & is.finite(r))
  if (length(good)) {
    aa = a[good]; bb = b[good]; rr = rep_len(r, n)[good]
    c = sqrt(1 - rr^2)
    lp = log_bvn_value(aa, bb, rr)
    ga = exp(stats::dnorm(aa, log = TRUE) +
      stats::pnorm((bb - rr*aa)/c, log.p = TRUE) - lp)
    gb = exp(stats::dnorm(bb, log = TRUE) +
      stats::pnorm((aa - rr*bb)/c, log.p = TRUE) - lp)
    gr = exp(-log(2*pi) - log(c) -
      (aa^2 - 2*rr*aa*bb + bb^2)/(2*c^2) - lp)
    logp[good] = lp
    score[good, ] = da[good, , drop = FALSE]*ga +
      db[good, , drop = FALSE]*gb + dr[good, , drop = FALSE]*gr
  }
  list(fit = exp(logp), gradient = score*exp(logp),
       log_fit = logp, log_gradient = score)
}

# 2. Prediction uncertainty

prediction_covariance = function(object, choice = NULL) {
  if (is.matrix(choice)) {
    covariance = choice
  } else {
    if (is.null(choice)) choice = object$control$vcov
    choice = match.arg(choice, c("observed", "opg", "sandwich"))
    covariance = stats::vcov(object, type = choice)
  }
  parameters = names(object$coefficients)
  if (!is.matrix(covariance) || !identical(rownames(covariance), parameters) ||
      !identical(colnames(covariance), parameters))
    stop("Prediction requires a named natural-scale covariance matrix", call. = FALSE)
  if (is.matrix(choice)) {
    if (!is.numeric(covariance) || is.complex(covariance) || any(!is.finite(covariance)))
      stop("Supplied covariance must contain finite numeric values", call. = FALSE)
    scale = max(abs(covariance), 0)
    if (scale > 0) {
      scaled = covariance / scale
      if (max(abs(scaled - t(scaled))) > 1e-10)
        stop("Supplied covariance must be symmetric", call. = FALSE)
      if (min(eigen((scaled + t(scaled))/2, symmetric = TRUE, only.values = TRUE)$values) < -1e-10)
        stop("Supplied covariance must be positive semidefinite", call. = FALSE)
    }
  }
  covariance
}

prediction_se = function(gradient, covariance) {
  if (all(is.finite(gradient)) && all(is.finite(covariance))) {
    variance = rowSums((gradient %*% covariance)*gradient)
    if (all(is.finite(variance)) && all(variance >= 0)) return(sqrt(variance))
  }
  # Unknown covariance entries involving a zero derivative do not invalidate
  # a target that does not depend on those parameters.
  apply(gradient, 1L, function(g) {
    if (anyNA(g)) return(NA_real_)
    ix = which(g != 0)
    if (!length(ix)) return(0)
    V = covariance[ix, ix, drop = FALSE]
    if (any(!is.finite(V))) return(NA_real_)
    variance = drop(g[ix] %*% V %*% g[ix])
    if (variance < -1e-10*max(1, sum(abs(V)))) return(NA_real_)
    sqrt(max(0, variance))
  })
}

# 3. Population and history integration

prediction_inputs = function(object, data, targets) {
  unused_design = function(design) list(matrix = matrix(0, nrow(data), length(design$columns),
    dimnames = list(NULL, design$columns)), offset = rep(0, nrow(data)))
  design = object$prepared$design
  px = if (all(targets %in% c('selection', 'selection_link'))) unused_design(design$outcome) else
    prediction_design(design$outcome, data)
  pw = if (all(targets %in% c('response', 'link'))) unused_design(design$selection) else
    prediction_design(design$selection, data)
  list(outcome = px, selection = pw)
}

population_predictions = function(object, data, targets, inputs = NULL) {
  if (is.null(inputs)) inputs = prediction_inputs(object, data, targets)
  model_predictions(object$coefficients, inputs$outcome$matrix, inputs$selection$matrix,
    inputs$outcome$offset, inputs$selection$offset, object$prepared$family,
    object$prepared$has_selection, targets = targets,
    structure = model_structure(object), model_id = object[["model_id"]] %||% object$prepared[["model_id"]])
}

prediction_result = function(object, data, type, conditioning, history) {
  inputs = prediction_inputs(object, data, type)
  px = inputs$outcome; pw = inputs$selection
  theta = object$coefficients
  components = population_predictions(object, data, type, inputs)
  result = components[[type]]
  if (conditioning == "population" || type %in% c("link", "selection_link")) return(result)
  if (is.null(history)) stop("History conditioning requires explicit history data", call. = FALSE)
  id = object$prepared$id
  if (!id %in% names(data) || anyNA(data[[id]]))
    stop("History prediction requires nonmissing target subject IDs", call. = FALSE)
  if (!is.data.frame(history)) history = as.data.frame(history)
  if (!nrow(history)) return(result)
  # History rows follow the fit's missing-data rule; padding has no meaning here.
  prepared = prepare_data(object$prepared$design$outcome$formula,
    if (object$prepared$has_selection) object$prepared$design$selection$formula else NULL,
    history, id, object$prepared$family,
    na_action = if (identical(object$control$na_action, "fail")) "fail" else "omit",
    design = object$prepared$design)
  rules = posterior_rule(theta, prepared, object$control, object$fixed)
  for (i in seq_len(nrow(data))) {
    j = match(data[[id]][i], prepared$id_values)
    if (is.na(j)) next
    rule = rules[[j]]; k = nrow(rule$nodes); p = length(theta)
    gradients = rule$node_gradient
    kernels = model_predictions(theta, px$matrix[rep(i, k), , drop = FALSE],
      pw$matrix[rep(i, k), , drop = FALSE], rep(px$offset[i], k), rep(pw$offset[i], k),
      object$prepared$family, object$prepared$has_selection,
      rule$nodes[, "u"], rule$nodes[, "v"],
      matrix(gradients[, 1L, ], k, p), matrix(gradients[, 2L, ], k, p),
      targets = c(type, if (type == 'selected') 'selection'),
      structure = model_structure(object), model_id = object[["model_id"]] %||% object$prepared[["model_id"]])
    target = kernels[[type]]
    log_weights = rule$log_weights
    score = rule$weight_score
    # Tilt the history posterior by the target selection probability. This
    # weighted selected mean equals E[D*Y* | history] / P(D=1 | history).
    if (type == "selected") {
      log_weights = log_weights + kernels$selection$log_fit
      score = score + kernels$selection$log_gradient
    }
    weights = exp(log_weights - max(log_weights))
    weights = weights/sum(weights)
    estimate = sum(weights*target$fit)
    # Differentiate a normalized weighted mean: E_w[target gradient] plus
    # E_w[(target - mean) * log-weight gradient], including moving nodes.
    gradient = colSums(target$gradient*weights) +
      colSums(score*(weights*(target$fit - estimate)))
    result$fit[i] = estimate; result$gradient[i, ] = gradient
  }
  result
}

# 4. Aggregation, row restoration and confidence bounds

finish_prediction = function(object, newdata, result, original = FALSE, se.fit = FALSE,
                              vcov = NULL, aggregate = FALSE, weights = NULL,
                              covariance = NULL, link = NULL, level = 0.95) {
  # A confidence interval (link not NULL) needs the standard error.
  standard_errors = se.fit || !is.null(link)
  names(result$fit) = rownames(newdata)
  rownames(result$gradient) = rownames(newdata)
  if (aggregate) {
    if (is.null(weights)) weights = rep(1, nrow(newdata))
    if (length(weights) != nrow(newdata) || any(!is.finite(weights)) ||
        any(weights < 0) || sum(weights) <= 0)
      stop('Averaging weights must be finite, nonnegative, and have positive sum', call. = FALSE)
    weights = weights/sum(weights)
    result$fit = sum(weights*result$fit)
    result$gradient = matrix(colSums(result$gradient*weights), 1L,
      dimnames = list(NULL, names(object$coefficients)))
  } else if (!is.null(weights)) stop('weights requires aggregate=TRUE', call. = FALSE)
  if (standard_errors) result$se.fit = prediction_se(result$gradient,
    if (is.null(covariance)) prediction_covariance(object, vcov) else covariance)
  if (original && !aggregate && inherits(object$prepared$na.action, 'exclude')) {
    rows = padded_rows(object$prepared, object$data)
    keep = match(object$prepared$keep, rows)
    restore = function(value) {
      if (is.matrix(value)) {
        ans = matrix(NA_real_, length(rows), ncol(value),
                     dimnames = list(NULL, colnames(value)))
        ans[keep, ] = value
      } else { ans = rep(NA_real_, length(rows)); ans[keep] = value }
      ans
    }
    result$fit = restore(result$fit); result$gradient = restore(result$gradient)
    names(result$fit) = rownames(object$data)[rows]
    rownames(result$gradient) = rownames(object$data)[rows]
    if (standard_errors) result$se.fit = restore(result$se.fit)
  }
  if (standard_errors) names(result$se.fit) = names(result$fit)
  if (!is.null(link)) result$fit = confidence_bounds(result$fit, result$se.fit, link, level)
  if (!se.fit) return(result$fit)
  result[c('fit', 'se.fit', 'gradient')]
}

# Delta-method confidence bounds on a link scale that keeps them in range:
# eta = g(m) has standard error |g'(m)| se, and the bounds are
# g^{-1}(eta -+ z se_eta). A target at the edge of its range (probability 0 or
# 1, mean 0) has no interior scale; both bounds equal the estimate there.
confidence_bounds = function(fit, se, link, level) {
  z = stats::qnorm((1 + level)/2)
  bounds = switch(link,
    identity = cbind(fit - z*se, fit + z*se),
    log = cbind(exp(log(fit) - z*se/fit), exp(log(fit) + z*se/fit)),
    probit = {
      eta = stats::qnorm(fit); half = z*se/stats::dnorm(eta)
      cbind(stats::pnorm(eta - half), stats::pnorm(eta + half))
    })
  edge = link != "identity" & !is.na(fit) & (fit <= 0 | (link == "probit" & fit >= 1))
  bounds[edge, ] = fit[edge]
  answer = cbind(fit = fit, lwr = bounds[, 1L], upr = bounds[, 2L])
  rownames(answer) = names(fit)
  answer
}

# Probabilities use the probit link, count means the log link, and Linear means
# and linear predictors the identity.
interval_link = function(type, family) {
  if (type %in% c("link", "selection_link")) return("identity")
  if (type == "selection" || family == "probit") return("probit")
  if (family %in% c("poisson", "pln")) "log" else "identity"
}

# 5. Public prediction and fitted-value methods

#' Predict outcome means and selection probabilities
#'
#' Compute a clearly defined mean, probability or linear predictor at the fitted
#' parameters. Population predictions average over person effects. History
#' predictions instead update each person's effect distribution using explicitly
#' supplied responses and selection indicators; they do not refit coefficients.
#'
#' @param object A fitted `panelselect` model from [panel_select()].
#' @param newdata Data frame containing the regressors, factor levels and offset
#'   variables needed by the target equations; response columns may be omitted.
#'   `NULL` (default) uses the retained fitting rows. Under `na_action = "exclude"`,
#'   default row predictions restore excluded rows as `NA`. History conditioning
#'   also requires the fitted ID column with nonmissing target IDs. New factor
#'   levels unsupported by the fitted design cause an error.
#' @param type One target: `"response"` (default potential-outcome mean),
#'   `"selection_probability"`, `"joint"`, `"selected_mean"`, or `"link"`.
#'   The table below distinguishes quantities that can differ substantially in
#'   selection models. Standalone models support `"response"` and `"link"` only.
#' @param equation Equation for `type = "link"`: `"outcome"` (default) or
#'   `"selection"`. Other target types already determine the equations they need;
#'   setting `equation = "selection"` for them is an error.
#' @param conditioning `"population"` (default) uses the prior person-effect
#'   distribution. `"history"` uses the posterior given the matching person's
#'   supplied history. History requires a panel model with GH/AGH integration;
#'   it is unavailable for GL, including PLN observation rules, and for links.
#' @param history Data frame of conditioning rows for `conditioning = "history"`.
#'   Supply the fitted ID column, selection indicator and selection regressors
#'   on all relevant rows, plus observed outcomes/outcome regressors on selected
#'   rows (all rows for standalone models). Include offset variables required by
#'   the formulas. `NULL` is an error for history prediction; no responses are
#'   taken implicitly from the target or training data. An empty history gives
#'   population predictions. With population conditioning, a supplied history
#'   is rejected rather than ignored.
#' @param se.fit Logical; `TRUE` returns parameter-estimation standard errors and
#'   natural-parameter gradients along with predictions. Default `FALSE` returns
#'   just predictions or the requested confidence-interval matrix.
#' @param vcov Covariance for SEs/intervals: `NULL` (default) uses the fit's choice;
#'   choose `"opg"`, `"observed"`, `"sandwich"`, or supply a named natural-scale
#'   covariance matrix. An observed/sandwich request may calculate information.
#'   This argument is used only when SEs or confidence intervals are requested.
#' @param aggregate Logical; default `FALSE` returns one prediction per row.
#'   `TRUE` averages the row predictions and their gradients, giving a single
#'   mean/probability and, if requested, its uncertainty.
#' @param weights Numeric averaging weights, one per row of `newdata`, used only
#'   with `aggregate = TRUE`. They must be finite, nonnegative and sum to a
#'   positive number; they are normalized internally. `NULL` gives equal row
#'   weights. These are prediction-averaging weights, not likelihood weights.
#' @param interval `"none"` (default) or `"confidence"` for confidence bounds on
#'   the estimated target. Bounds include parameter uncertainty only;
#'   `"prediction"` intervals for realized outcomes are not implemented.
#' @param level Confidence level strictly between zero and one when
#'   `interval = "confidence"`; default `0.95`.
#' @param ... `predict()` rejects unused arguments. `fitted()` forwards its
#'   arguments to `predict()`, so the same targets and options are available.
#'
#' @section Choose the target:
#' Write `D` for selection and `Y*` for the potential outcome before selection.
#' Means below condition on the supplied regressors and any explicit history.
#'
#' | `type` | Quantity | Use |
#' |---|---|---|
#' | `"response"` | `E[Y*]` | Mean outcome if every target row could be observed; this is also the default for `fitted()`. |
#' | `"selection_probability"` | `P(D = 1)` | Probability that the outcome is observed or the person participates. |
#' | `"joint"` | `E[D * Y*]` | Mean observed contribution when nonparticipation contributes zero. |
#' | `"selected_mean"` | `E[Y* | D = 1]` | Mean among selected/participating targets; equals the joint mean divided by selection probability. |
#' | `"link"` | Fixed linear predictor including offsets | Inspect the chosen equation's regression index, without integrating person effects or applying the response link. |
#'
#' Binary-outcome response means are probabilities. Linear-outcome means use
#' outcome units, while Poisson/PLN means are expected counts, including exposure
#' offsets when specified. Each target represents a fresh within-person shock.
#'
#' @section Choose the conditioning information:
#' - **Population:** integrate over the model's person-effect distribution. This
#'   is appropriate for a new person or a population-average target. It is
#'   available after GH, AGH and GL fitting and uses analytic target formulas.
#' - **History:** integrate over each person's posterior given exactly the rows
#'   supplied in `history`. IDs absent from usable history receive population
#'   predictions. The current GH/AGH history calculation uses the fit's stored
#'   integration settings; there is no prediction-specific method override.
#'
#' The function does not infer chronology, past-only windows or leave-one-out
#' samples. For forecasting, supply only information available before the target.
#' Including a target's own observed response informs its person effects and
#' produces a fresh-wave replicated prediction, not a held-out forecast or a
#' prediction of that same realized shock. An absent history row is not an
#' observed `D = 0` event. Missing required history values follow the fitted
#' `na_action`: `"fail"` stops, otherwise incomplete history rows are omitted.
#'
#' For selected means the target selection event also reweights the history
#' posterior. The result is a ratio of integrated moments, not an unweighted
#' average of node-specific selected means.
#'
#' @section Uncertainty and aggregation:
#' - SEs use the delta method with uncertainty in the fitted parameters. For
#'   history targets, derivatives include posterior weights and moving nodes.
#' - Confidence bounds use a probit scale for probabilities, a log scale for
#'   count means, and the identity scale for Linear means and linear predictors,
#'   then transform back. Aggregation occurs before this interval calculation.
#' - A supplied matrix must name all natural parameters, be finite and symmetric,
#'   and be positive semidefinite (singular matrices are allowed). Symmetry and
#'   eigenvalue checks allow relative roundoff of `1e-10`.
#' - Unavailable covariance produces unavailable SEs; a different estimator is
#'   never substituted. These intervals do not include realized-outcome variation.
#'
#' @return Normally a named numeric vector; with `aggregate = TRUE`, a single
#'   estimate. `interval = "confidence"` returns a matrix with `fit`, `lwr` and
#'   `upr`. With `se.fit = TRUE`, the result is a list containing that `fit`
#'   vector/matrix, a numeric `se.fit`, and a `gradient` matrix with natural
#'   parameters in columns. Fixed-parameter gradient columns do not add variance.
#' @section Fitted values versus residuals:
#' `fitted()` defaults to the potential-outcome mean, whereas [residuals()] uses
#' the selected mean by default for selection models. Therefore adding their
#' default results does not generally reconstruct observed selected outcomes.
#' @family prediction functions
#' @seealso [random_effects()], [residuals.panelselect()], [panelselect-inference]
#' @examples
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#' fit = panel_select(y ~ x, observed ~ x + z, dat, "person", "linear")
#'
#' targets = data.frame(person = 1:2, x = 0, z = 0)
#' predict(fit, targets, type = "selection_probability", interval = "confidence")
#' # Condition on each person's observed rows.
#' predict(fit, targets, type = "selected_mean", conditioning = "history",
#'         history = dat[dat$person %in% 1:2, ])
#' predict(fit, aggregate = TRUE, se.fit = TRUE)$fit
#' @export
predict.panelselect = function(object, newdata = NULL,
    type = c("response", "selection_probability", "joint", "selected_mean", "link"),
    equation = c("outcome", "selection"),
    conditioning = c("population", "history"), history = NULL, se.fit = FALSE,
    vcov = NULL, aggregate = FALSE, weights = NULL, interval = c("none", "confidence"),
    level = 0.95, ...) {
  reject_unused(list(...), "predict", "see ?predict.panelselect for its arguments")
  if (identical(interval, "prediction"))
    stop('interval = "prediction" is not available: uncertainty concerns the estimated mean or ',
         'probability, not a realized outcome; use interval = "confidence"', call. = FALSE)
  interval = match.arg(interval)
  if (interval == "confidence" && (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1))
    stop("level must lie strictly between zero and one", call. = FALSE)
  object = current_fit(object)
  earlier = c(selection = "selection_probability", selected = "selected_mean")
  if (is.character(type) && length(type) == 1L && type %in% names(earlier))
    stop('type = "', type, '" is now "', earlier[[type]], '"', call. = FALSE)
  type = match.arg(type); conditioning = match.arg(conditioning)
  # The prediction components keep their short internal names.
  type = switch(type, selection_probability = "selection", selected_mean = "selected", type)
  if (uses_nonhistory_quadrature(object$control, object$prepared) &&
      (conditioning == "history" || !is.null(history)))
    stop("History prediction is unavailable for GL integration, including PLN observation rules", call. = FALSE)
  if (model_structure(object) == "cross_sectional" &&
      (conditioning == "history" || !is.null(history)))
    stop("History conditioning is not defined for cross-sectional models", call. = FALSE)
  equation = match.arg(equation)
  if (equation == "selection") {
    if (type != "link") stop("equation applies only to type='link'", call. = FALSE)
    if (!object$prepared$has_selection) stop("The model has no selection equation", call. = FALSE)
    type = "selection_link"
  }
  if (!object$prepared$has_selection && type %in% c("selection", "joint", "selected"))
    stop("This prediction type requires a selection equation", call. = FALSE)
  # A supplied history is never ignored silently; link predictors have no
  # person-effect integral to condition.
  if (conditioning == "population" && !is.null(history))
    stop("history is used only with conditioning = 'history'", call. = FALSE)
  if (conditioning == "history" && type %in% c("link", "selection_link"))
    stop("Link predictions do not depend on history; use conditioning = 'population'", call. = FALSE)
  original = is.null(newdata)
  if (original) newdata = object$data[object$prepared$keep, , drop = FALSE]
  newdata = as.data.frame(newdata)
  result = prediction_result(object, newdata, type, conditioning, history)
  finish_prediction(object, newdata, result, original, se.fit, vcov, aggregate, weights,
    link = if (interval == "confidence") interval_link(type, object$prepared$family), level = level)
}

#' @rdname predict.panelselect
#' @export
fitted.panelselect = function(object, ...) stats::predict(object, ...)

# 6. Residuals on the fitting sample

#' Calculate response residuals for observed outcomes
#'
#' Subtract a predicted mean from each observed outcome on the fitting sample.
#' The default target is the selected mean in a selection model and the response
#' mean in a standalone model. Unselected outcomes have no residual.
#'
#' @param object A fitted `panelselect` model. Residuals use its retained fitting
#'   rows; this method is not a new-data prediction interface.
#' @param ... Options forwarded to [predict.panelselect()]. `type` selects the
#'   mean to subtract, not a residual standardization: default `"selected_mean"`
#'   for selection models and `"response"` otherwise; the other supported
#'   prediction targets retain their usual meanings. `conditioning = "history"`
#'   requires explicit `history` and supported GH/AGH rules. Requests for new
#'   data, aggregated output, SE output or interval matrices are not supported.
#'   There are no Pearson or deviance residuals.
#'
#' @section Common choices:
#' - `residuals(fit)` uses the population mean appropriate to the observed
#'   outcomes: selected mean for selection models, response mean otherwise.
#' - `residuals(fit, type = "response")` uses the potential-outcome mean;
#'   selection can make these residuals systematically nonzero among participants.
#' - `residuals(fit, conditioning = "history", history = dat)` uses person effects
#'   informed by `dat`. If it includes each row's own outcome, these are
#'   descriptive residuals and should not be interpreted as forecast errors.
#'
#' `fitted()` has a different default target in selection models. To reconstruct
#' selected outcomes, use `predict(fit, type = "selected_mean") + residuals(fit)`.
#'
#' @return A numeric vector in fitting-row order, with `NA` for unselected rows.
#'   Under `na_action = "exclude"`, excluded fitting rows are also restored as
#'   `NA`; `"omit"` does not restore dropped rows.
#' @family prediction functions
#' @seealso [predict.panelselect()], [random_effects()]
#' @examples
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#' fit = panel_select(y ~ x, observed ~ x + z, dat, "person", "linear")
#' summary(residuals(fit))
#' # Residuals from each person's posterior mean given their own rows.
#' head(residuals(fit, conditioning = "history", history = dat))
#' @export
residuals.panelselect = function(object, ...) {
  object = current_fit(object)
  args = list(...)
  prediction_call = as.call(c(list(quote(predict.panelselect), object = quote(object)), args))
  matched = match.call(predict.panelselect, prediction_call)
  if (!is.null(matched$newdata))
    stop("Residuals use the retained fitting rows; newdata is not supported", call. = FALSE)
  if (is.null(matched$type)) args$type = if (object$prepared$has_selection) "selected_mean" else "response"
  # The matched type covers positional and abbreviated arguments.
  type = matched$type %||% args$type
  targets = eval(formals(predict.panelselect)$type)
  if (!is.character(type) || length(type) != 1L || is.na(pmatch(type, targets)))
    stop("residuals() returns response residuals, y minus a predicted mean; type selects that mean (",
         paste0('"', targets, '"', collapse = ", "), "), not a residual type such as \"pearson\"",
         call. = FALSE)
  fitted = do.call(stats::predict, c(list(object = object), args))
  exclude = inherits(object$prepared$na.action, "exclude")
  rows = padded_rows(object$prepared, object$data)
  expected = if (exclude) length(rows) else length(object$prepared$keep)
  if (!is.numeric(fitted) || length(fitted) != expected)
    stop("Residuals require unaggregated predictions without se.fit", call. = FALSE)
  y = object$prepared$y[object$prepared$inverse_index]
  d = object$prepared$d[object$prepared$inverse_index]
  y[d == 0] = NA_real_
  if (exclude) {
    full = rep(NA_real_, length(rows)); full[match(object$prepared$keep, rows)] = y; y = full
  }
  y - fitted
}

# Rows that na_action = "exclude" pads: every data row, or the rows of a subset.
padded_rows = function(prepared, data) {
  subset = prepared$subset
  if (is.null(subset)) seq_len(nrow(data)) else which(is.na(subset) | subset)
}


# 7. Posterior person-effect summaries

#' Summarize person effects conditional on observed history
#'
#' Calculate posterior means and SDs of person-specific random-intercept
#' contributions at the fitted parameters. These summaries describe conditional
#' heterogeneity in each equation, rather than re-estimating its coefficients.
#'
#' @param object A fitted panel `panelselect` model using GH/AGH integration.
#'   Cross-sectional models have no person effects. GL at the person-effect or
#'   PLN observation level is not supported by this function.
#' @param history Data frame of conditioning rows, including the fitted ID and
#'   the response, selection, regressor and offset columns required by the
#'   original formulas. Required outcome data apply only to selected rows in
#'   selection models. `NULL` (default) uses the retained fitting rows, unlike
#'   [predict()] history conditioning, which always requires explicit history.
#'   Supply a restricted history for an earlier information set. Missing
#'   required values stop under `na_action = "fail"` and are omitted otherwise.
#'
#' @section Interpretation and available summaries:
#' | Column | Meaning |
#' |---|---|
#' | Fitted ID-column name | Person whose supplied usable rows define the posterior. |
#' | `selection_mean`, `selection_sd` | Posterior mean and SD of the selection random intercept; included only for selection models. |
#' | `outcome_mean`, `outcome_sd` | Posterior mean and SD of the outcome random intercept. |
#'
#' The intercepts include their fitted person-effect SD multipliers. For count
#' models they contribute to the log rate, not directly to the expected count.
#' The posterior SDs measure remaining uncertainty about a person's effect at
#' fixed fitted parameters; they are not parameter-estimation SEs or confidence
#' intervals for regression coefficients.
#'
#' Exactly the supplied history is used, without automatic date filtering.
#' Default summaries use the person's own fitting responses and are descriptive.
#' Only IDs represented in usable history are returned; an empty history returns
#' an empty data frame with the appropriate columns. For predictions of people
#' without observed history, use [predict()] and its population fallback.
#'
#' @return A data frame with one row per history person and the columns above.
#'   A structurally fixed-zero effect has zero posterior mean and SD.
#' @family prediction functions
#' @seealso [predict.panelselect()], [panel_select()]
#' @examples
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#' fit = panel_select(y ~ x, observed ~ x + z, dat, "person", "linear")
#' head(random_effects(fit))
#' @export
random_effects = function(object, history = NULL) {
  object = current_fit(object)
  if (uses_nonhistory_quadrature(object$control, object$prepared))
    stop("Random effects are unavailable for GL integration, including PLN observation rules", call. = FALSE)
  if (model_structure(object) == "cross_sectional")
    stop("Cross-sectional models have no person random effects", call. = FALSE)
  if (is.null(history)) history = object$data[object$prepared$keep, , drop = FALSE]
  history = as.data.frame(history)
  id = object$prepared$id
  if (!id %in% names(history)) stop("History must contain the subject ID column", call. = FALSE)
  columns = c(if (object$prepared$has_selection) c("selection_mean", "selection_sd"),
               "outcome_mean", "outcome_sd")
  if (!nrow(history)) {
    answer = history[, id, drop = FALSE]
    for (name in columns) answer[[name]] = numeric()
    return(answer)
  }
  prepared = prepare_data(object$prepared$design$outcome$formula,
    if (object$prepared$has_selection) object$prepared$design$selection$formula else NULL,
    history, id, object$prepared$family,
    if (identical(object$control$na_action, "fail")) "fail" else "omit",
    design = object$prepared$design)
  rules = posterior_rule(object$coefficients, prepared, object$control, object$fixed)
  answer = data.frame(prepared$id_values); names(answer) = id
  for (name in columns) answer[[name]] = NA_real_
  for (i in seq_along(rules)) {
    weights = exp(rules[[i]]$log_weights)
    weights = weights/sum(weights)
    for (equation in c(if (prepared$has_selection) "selection", "outcome")) {
      latent = if (equation == "selection") "u" else "v"
      scale = object$coefficients[[paste0("sd_re_", equation)]]
      values = scale*rules[[i]]$nodes[, latent]
      mean = sum(weights*values)
      answer[[paste0(equation, "_mean")]][i] = mean
      answer[[paste0(equation, "_sd")]][i] = sqrt(sum(weights*(values-mean)^2))
    }
  }
  answer
}
