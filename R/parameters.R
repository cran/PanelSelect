# Natural parameters and their optimizer coordinates.
# parameter_layout() names the natural parameters: selection and
# outcome coefficients, then the model's SDs and correlations. The optimizer
# works with log SDs and atanh correlations (parameter_to_optimizer(),
# parameter_from_optimizer(), with the Jacobian parameter_jacobian()).
# resolve_fixed() validates restrictions, and initial_parameters() builds
# pooled or separate random-effects starts (see panel_select()).
# 1. Natural parameter names and restrictions

parameter_layout = function(prepared) {
  model = resolve_model(prepared$family, prepared$has_selection,
    prepared$structure, prepared$model_id)
  names = c(if (prepared$has_selection && ncol(prepared$w))
               paste0("selection.", colnames(prepared$w)),
             if (ncol(prepared$x)) paste0("outcome.", colnames(prepared$x)),
             model$parameters)
  data.frame(name = names,
             type = ifelse(startsWith(names, "sd_"), "sd",
                           ifelse(startsWith(names, "cor_"), "correlation", "coefficient")),
             fixed = FALSE, stringsAsFactors = FALSE)
}

# Earlier person-effect parameter names. New inputs report their replacement;
# current_fit() translates saved fits without changing data identifiers.
renamed_parameters = c(sd_selection = "sd_re_selection", sd_outcome = "sd_re_outcome",
  sd_id_selection = "sd_re_selection", sd_id_outcome = "sd_re_outcome", cor_id = "cor_re")

current_names = function(x) {
  earlier = x %in% names(renamed_parameters)
  x[earlier] = renamed_parameters[x[earlier]]
  x
}

validate_named_parameters = function(values, info, label, allow_zero = FALSE) {
  if (is.null(values)) return(numeric())
  if (!is.numeric(values) || is.complex(values) || !is.null(dim(values)) ||
      (length(values) && (is.null(names(values)) || anyNA(names(values)) ||
                         any(names(values) == "") || anyDuplicated(names(values)))) ||
      any(!is.finite(values)))
    stop(label, " must be a finite named numeric vector", call. = FALSE)
  earlier = intersect(names(values), names(renamed_parameters))
  if (length(earlier))
    stop(paste0(earlier, " is now ", renamed_parameters[earlier], collapse = "; "), " (in ", label, ")",
         call. = FALSE)
  if (any(!names(values) %in% info$name))
    stop("Unknown ", label, " parameter(s): ",
         paste(setdiff(names(values), info$name), collapse = ", "), "; valid names: ",
         paste(info$name, collapse = ", "), call. = FALSE)
  type = info$type[match(names(values), info$name)]
  if (any(type == "sd" & (values < 0 | (!allow_zero & values == 0))))
    stop(label, " standard deviations must be ", if (allow_zero) "nonnegative" else "positive",
         call. = FALSE)
  if (any(type == "correlation" & abs(values) >= 1))
    stop(label, " correlations must lie strictly between -1 and 1", call. = FALSE)
  values
}

resolve_fixed = function(fixed, info, prepared) {
  fixed = validate_named_parameters(fixed, info, "fixed", allow_zero = TRUE)
  model = resolve_model(prepared$family, prepared$has_selection,
    prepared$structure, prepared$model_id)
  if (any(fixed[intersect(names(fixed), model$positive_scales)] == 0))
    stop("The residual standard deviation must be positive", call. = FALSE)
  absent = names(fixed)[fixed == 0]
  redundant = names(Filter(function(scales) any(scales %in% absent),
                            model$zero_correlations))
  for (name in redundant) {
    if (name %in% names(fixed) && fixed[[name]] != 0)
      stop(name, " is unidentified at the fixed zero scale (",
           paste0(intersect(model$zero_correlations[[name]], absent), " = 0", collapse = ", "),
           ") and must equal zero; it is fixed at zero automatically when omitted", call. = FALSE)
    fixed[name] = 0
  }
  fixed
}

# 2. Optimizer coordinates and their Jacobian

# The transforms below index by position in parameter_layout order.
check_parameter_order = function(theta, info)
  if (!identical(names(theta), info$name))
    stop("Parameter values do not follow the model's parameter order", call. = FALSE)

parameter_to_optimizer = function(theta, info) {
  check_parameter_order(theta, info)
  free = !info$fixed
  value = theta[free]
  type = info$type[free]
  value[type == "sd"] = log(value[type == "sd"])
  value[type == "correlation"] = atanh(value[type == "correlation"])
  value
}

parameter_from_optimizer = function(eta, template, info, check = TRUE) {
  if (check) check_parameter_order(template, info)
  free = !info$fixed
  type = info$type[free]
  value = eta
  value[type == "sd"] = exp(value[type == "sd"])
  value[type == "correlation"] = tanh(value[type == "correlation"])
  template[free] = value
  template
}

parameter_jacobian = function(theta, info, check = TRUE) {
  if (check) check_parameter_order(theta, info)
  value = stats::setNames(rep(1, nrow(info)), info$name)
  value[info$type == "sd"] = theta[info$type == "sd"]
  value[info$type == "correlation"] = 1 - theta[info$type == "correlation"]^2
  value[info$fixed] = 0
  value
}

# 3. Explicit starts, restrictions and initialization precedence

validate_initialization = function(init) {
  choices = c(default = "default", re = "default", pooled = "pooled",
    cross_sectional = "pooled", zero = "zero", unif = "unif", norm = "norm")
  if (!is.character(init) || length(init) != 1L || !is.null(dim(init)) ||
      is.na(init) || !init %in% names(choices))
    stop("Unsupported initialization method: init must be one of ",
      paste0('"', names(choices), '"', collapse = ", "), call. = FALSE)
  unname(choices[init])
}

initial_parameters = function(prepared, info, start, fixed, control) {
  start = validate_named_parameters(start, info, "start", allow_zero = TRUE)
  both = intersect(names(start), names(fixed))
  conflict = both[start[both] != fixed[both]]
  if (length(conflict))
    stop("A supplied start conflicts with a fixed parameter: ", paste(conflict, collapse = ", "),
         call. = FALSE)
  if (any(start[names(start) %in% info$name[info$type == "sd" & !info$fixed]] <= 0))
    stop("Starting free standard deviations must be positive", call. = FALSE)
  theta = stats::setNames(rep(0, nrow(info)), info$name)
  theta[info$type == "sd"] = 1
  theta[names(fixed)] = fixed
  theta[names(start)] = start
  free_names = info$name[!info$fixed]
  scale = function(name) if (name %in% names(theta)) theta[[name]] else 0
  if (all(free_names %in% names(start))) return(theta)

  if (control$init %in% c("zero", "unif", "norm")) {
    eta = stats::setNames(rep(0, length(free_names)), free_names)
    if (control$init == "unif") eta[] = stats::runif(length(eta))
    if (control$init == "norm") eta[] = stats::rnorm(length(eta))
    theta = parameter_from_optimizer(eta, theta, info)
  } else {
    pooled = function(x, y, offset, family, scale = 0) {
      if (!ncol(x)) return(numeric())
      fit = if (family == "gaussian") stats::lm.fit(x, y - offset) else {
        link = if (family == "probit") stats::binomial("probit") else stats::poisson()
        multiplier = if (family == "probit") sqrt(1 + scale^2) else 1
        offset = if (family == "probit") offset / multiplier else offset
        value = suppressWarnings(stats::glm.fit(x, y, family = link, offset = offset))
        value$coefficients = value$coefficients * multiplier
        value
      }
      value = fit$coefficients
      value[!is.finite(value)] = 0
      value
    }
    if (prepared$has_selection && ncol(prepared$w)) {
      theta[paste0("selection.", colnames(prepared$w))] =
        pooled(prepared$w, prepared$d, prepared$offset_selection, "probit", scale("sd_re_selection"))
    }
    take = prepared$selected
    if (any(take) && ncol(prepared$x)) theta[paste0("outcome.", colnames(prepared$x))] =
      pooled(prepared$x[take, , drop = FALSE], prepared$y[take],
             prepared$offset_outcome[take], prepared$family, scale("sd_re_outcome"))
    # Linear scales start from the pooled residual SD of the selected outcomes,
    # its variance split equally between the free scales not supplied as starts.
    linear = setdiff(intersect(c("sd_re_outcome", "sd_error"), free_names), names(start))
    if (prepared$family == "gaussian" && any(take) && length(linear)) {
      residual = prepared$y[take] - prepared$offset_outcome[take]
      rank = 0L
      if (ncol(prepared$x)) {
        fit = stats::lm.fit(prepared$x[take, , drop = FALSE], residual)
        residual = fit$residuals; rank = fit$rank
      }
      s = sqrt(sum(residual^2)/max(1, length(residual) - rank))
      if (is.finite(s) && s > 0) theta[linear] = s/sqrt(length(linear))
    }
    if (control$init %in% c("default", "re"))
      theta = re_initial_parameters(prepared, theta, start, fixed, control)
  }
  theta[names(fixed)] = fixed
  theta[names(start)] = start
  if (any(!is.finite(theta))) stop("Automatic starts are not finite; supply start", call. = FALSE)
  theta
}

# 4. Separate random-effects auxiliary fits

# Probit RE seeds selection; the matching RE family seeds the selected outcome.
# Reuse the retained design matrices: rebuilding formulas on a selected subset
# could change factor columns or data-dependent polynomial/spline bases.
re_initial_parameters = function(prepared, theta, start, fixed, control) {
  if (!identical(prepared$structure, "panel") || !prepared$has_selection) return(theta)
  supplied = union(names(start), names(fixed))
  initial = theta; initial[names(fixed)] = fixed; initial[names(start)] = start
  for (equation in c("selection", "outcome")) {
    selection = equation == "selection"
    family = if (selection) "probit" else prepared$family
    design = if (selection) prepared$w else prepared$x
    scale = if (selection) "sd_re_selection" else "sd_re_outcome"
    scales = c(scale, if (family %in% c("gaussian", "pln")) "sd_error")
    joint_names = c(if (ncol(design)) paste0(equation, ".", colnames(design)), scales)
    if (all(joint_names %in% supplied)) next
    rows = if (selection) seq_len(prepared$n_rows) else which(prepared$selected)
    # Fixed-zero probit/Poisson effects reduce to the pooled fit. Without any
    # repeated usable rows, probit cannot identify its free person SD; Gaussian
    # and PLN cannot separate two free SDs. Keep pooled starts in these cases.
    zero = function(name) name %in% names(fixed) && fixed[[name]] == 0
    if (family %in% c("probit", "poisson") && zero(scale)) next
    if (family == "pln" && zero(scale) && zero("sd_error")) next
    repeated = any(tabulate(prepared$group[rows], prepared$n_panels) > 1L)
    if (!repeated && ((family == "probit" && !scale %in% names(fixed)) ||
        (family %in% c("gaussian", "pln") && !any(scales %in% names(fixed))))) next
    # Quadrature follows the equation; joint optimizer options do not control
    # these auxiliary fits. A PLN outcome also keeps its observation rule.
    # Fresh controls deliberately carry no retained mesh or optimizer state.
    auxiliary_control = panel_control(quadrature = control$quadrature,
      nodes = rep_len(control$nodes, 2L)[if (selection) 1L else 2L],
      mode_maxit = control$mode_maxit, mode_tol = control$mode_tol,
      optimizer = "optim", method = "BFGS", quadrature_options = control$quadrature_options,
      observation_quadrature = control$observation_quadrature,
      observation_nodes = control$observation_nodes, observation_options = control$observation_options)
    automatic = tryCatch({
      aux_fixed = fixed[intersect(joint_names, names(fixed))]
      if (family == "gaussian") {
        aux_names = joint_names
        fit = linear_re_initial_fit(prepared, initial[joint_names], aux_fixed, auxiliary_control)
      } else {
        data = data.frame(.response = if (selection) prepared$d[rows] else prepared$y[rows],
          .id = prepared$group[rows],
          .offset = if (selection) prepared$offset_selection[rows] else prepared$offset_outcome[rows])
        formula = .response ~ 0 + offset(.offset)
        if (ncol(design)) {
          data$.design = I(design[rows, , drop = FALSE])
          formula = .response ~ 0 + .design + offset(.offset)
        }
        columns = colnames(stats::model.matrix(formula, data))
        aux_names = c(if (length(columns)) paste0("outcome.", columns), "sd_re_outcome",
          if (family == "pln") "sd_error")
        names(aux_fixed) = aux_names[match(names(aux_fixed), joint_names)]
        # Auxiliary covariance/score screens do not describe the joint fit.
        # Accept starts by backend status; the joint fit reports its own screens.
        fit = withCallingHandlers(panel_select(formula, NULL, data, ".id", family,
          start = stats::setNames(initial[joint_names], aux_names),
          fixed = aux_fixed, control = auxiliary_control, init = "pooled", vcov = "opg"),
          warning = function(warning) {
            message = conditionMessage(warning)
            if (startsWith(message, "The requested covariance estimate is unavailable;") ||
                startsWith(message, "Optimizer did not report convergence:") ||
                startsWith(message, "Optimizer stopped with Q_OPG = "))
              invokeRestart("muffleWarning")
          })
      }
      if (!isTRUE(fit$optimizer$converged)) stop("optimizer did not converge")
      value = fit$coefficients
      if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
          length(value) != length(aux_names) || anyDuplicated(names(value)) ||
          !setequal(names(value), aux_names) || any(!is.finite(value)))
        stop("auxiliary coefficients are missing, nonfinite, or incorrectly named")
      if (length(aux_fixed) && any(value[names(aux_fixed)] != aux_fixed))
        stop("auxiliary estimates do not preserve fixed parameters")
      value = stats::setNames(value[aux_names], joint_names)
      free_scales = setdiff(scales, names(fixed))
      if (any(value[scales] < 0) || any(value[free_scales] <= 0) ||
          (family == "gaussian" && value[["sd_error"]] <= 0))
        stop("auxiliary standard deviations are invalid")
      if (!is.null(fit$loglik) && (length(fit$loglik) != 1L || !is.finite(fit$loglik)))
        stop("auxiliary log likelihood is not finite")
      value
    }, error = function(error) {
      warning(equation, " RE initialization failed; using pooled starts: ",
        conditionMessage(error), call. = FALSE)
      NULL
    })
    if (!is.null(automatic)) theta[names(automatic)] = automatic
  }
  theta
}
