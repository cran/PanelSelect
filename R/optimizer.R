# Optimizer coordinates and backends: the evaluator that caches one likelihood
# per parameter point, BFGS scaling, the optim/maxLik call and status codes.

# 1. Likelihood evaluator and information endpoint ---------------------------

fit_evaluator = function(template, info, prepared, control, fixed) {
  control = current_control(control)
  stop_removed_settings(control)
  # Evaluations keep the template's parameter order, checked here once.
  check_parameter_order(template, info)
  information_scales = NULL
  staged = uses_staged_integration(prepared, control, fixed)
  moving_cut = uses_gl_quadrature(control) && !staged
  if (moving_cut) {
    prepared = prepare_model_data(prepared)
    prepared$split_cache = prepare_split_panels(prepared)
    prepared$split_workspace = new_split_workspace(prepared)
  }
  if (staged) {
    prepared = prepare_model_data(prepared)
    # Evaluator-local reuse of immutable mesh work; never stored with the fit.
    prepared$staged_workspace = new.env(parent = emptyenv())
  }
  if (!moving_cut && (!staged || control$quadrature %in% c('gh', 'agh'))) {
    prepared = prepare_integration(prepared, control)
    control = prepared$integration_cache$control
  }
  modes = NULL
  if (!moving_cut && control$quadrature == 'agh') {
    reuse = prepare_mode_reuse(prepared)
    prepared = reuse$prepared; modes = reuse$state
  }
  cache = new.env(parent = emptyenv())
  cache$eta = NULL
  cache$value = NULL
  cache$rejected = 0L
  cache$evaluations = 0L
  cache$best_valid = NULL
  cache$rejections = list(total = 0L, by_class = integer(), by_reason = integer(), last = NULL)
  free = !info$fixed
  evaluate = function(eta, trial = TRUE) {
    if (identical(eta, cache$eta)) {
      if (!trial && !cache$value$valid) stop(cache$value$error)
      return(cache$value)
    }
    if (!is.null(modes)) modes$pending = list()
    # Distinct likelihood evaluations, including rejected trials.
    cache$evaluations = cache$evaluations + 1L
    value = tryCatch({
      theta = parameter_from_optimizer(eta, template, info, check = FALSE)
      if (any(!is.finite(theta)) || any(theta[info$type == "sd" & free] <= 0) ||
          any(abs(theta[info$type == "correlation"]) >= 1))
        integration_error("Optimizer trial parameters are not representable")
      result = evaluate_model(theta, prepared, control, fixed, score = TRUE)
      if (length(result$loglik) != prepared$n_panels || any(!is.finite(result$loglik)))
        integration_error("The panel likelihood is not finite")
      if (!is.matrix(result$score) ||
          !identical(dim(result$score), c(prepared$n_panels, nrow(info))) ||
          any(!is.finite(result$score[, free, drop = FALSE])))
        integration_error("The panel score is unavailable or has invalid dimensions")
      colnames(result$score) = info$name
      result$score[, !free] = 0
      jacobian = parameter_jacobian(theta, info, check = FALSE)[free]
      result$optimizer_score = sweep(result$score[, free, drop = FALSE], 2L, jacobian, `*`)
      if (moving_cut || staged) result$reference_optimizer_score = sweep(result$reference_score[, free, drop = FALSE],
        2L, jacobian, `*`)
      result$theta = theta
      result$valid = TRUE
      result
    }, panelselect_integration_error = function(error) {
      cache$rejected = cache$rejected + 1L
      # Keep compact rejection evidence, not every trial or its model arrays.
      kind = class(error)[1L]; reason = conditionMessage(error)
      keys = c(by_class = kind, by_reason = reason)
      for (field in names(keys)) {
        key = keys[[field]]
        if (!key %in% names(cache$rejections[[field]])) cache$rejections[[field]][key] = 0L
        cache$rejections[[field]][key] = cache$rejections[[field]][key] + 1L
      }
      cache$rejections$total = cache$rejected
      cache$rejections$last = list(class = kind, message = reason, diagnostics = error[["diagnostics"]])
      if (inherits(error, "panelselect_mesh_required")) stop(error)
      list(valid = FALSE, error = error,
           loglik = rep(-Inf, prepared$n_panels),
           optimizer_score = matrix(0, prepared$n_panels, sum(free)))
    })
    cache$eta = eta
    cache$value = value
    # A staged restart starts from the best valid point of its stage. One stage
    # has one immutable mesh, so its likelihood values are comparable.
    if (staged && value$valid && (is.null(cache$best_valid) ||
        sum(value$loglik) > cache$best_valid$loglik))
      cache$best_valid = list(eta = eta, theta = value$theta, loglik = sum(value$loglik))
    if (!is.null(modes)) {
      if (value$valid && is.null(modes$anchor))
        modes$seeds[names(modes$pending)] = modes$pending
      modes$pending = list()
    }
    if (!trial && !value$valid) stop(value$error)
    value
  }
  # Information differences use one fitted-point seed for both signs and every
  # parameter. Restore it afterwards, even when differentiation fails.
  anchor = function(eta) {
    if (is.null(modes)) return(function() invisible(NULL))
    evaluate(eta, trial = FALSE)
    previous = modes$anchor; seeds = modes$seeds
    modes$anchor = seeds
    function() {
      modes$seeds = seeds; modes$anchor = previous; modes$pending = list()
      cache$eta = cache$value = NULL
      invisible(NULL)
    }
  }
  answer = list(evaluate = evaluate, rejected = function() cache$rejected,
       rejections = function() cache$rejections, evaluations = function() cache$evaluations,
       anchor = anchor,
       best_valid = function() cache$best_valid, staged = staged,
       # Only immutable design scales are cached. The Gaussian residual SD is
       # taken from this endpoint, never from the starting template or an earlier
       # information request. OPG-only fits do not calculate these steps.
       information_steps = function(eta) {
         if (is.null(information_scales)) information_scales <<- design_scales(prepared)
         theta = parameter_from_optimizer(eta, template, info, check = FALSE)
         information_steps(info, prepared, theta, information_scales)
       },
       mode_counts = function() if (!is.null(modes)) modes$counts else NULL)
  if (moving_cut) {
    answer$information_gate_tolerance = control$quadrature_options$information_tol
    answer$information_gate_basis = paste("Analytical finite moving-cut GL scores including moving nodes and interval weights;",
      "native observed information uses canonical/predecessor orders and two score-difference steps.")
  }
  if (staged) {
    answer$information_gate_tolerance = staged_information_tolerance(control, prepared$family)
    answer$information_gate_basis = paste("Analytical scores of immutable positive meshes;",
      "canonical/reference rules and two score-difference steps share each retained stage.")
  }
  answer
}

# 2. Backend calls and convergence status ------------------------------------

optimize_model = function(eta, evaluator, control, options = optimizer_options(control)) {
  objective = function(value) -sum(evaluator$evaluate(value)$loglik)
  gradient = function(value) -colSums(evaluator$evaluate(value)$optimizer_score)
  if (!length(eta)) {
    result = list(par = eta, counts = c("function" = 1L, gradient = 1L),
                   convergence = 0L, message = "All parameters fixed")
  } else if (control$optimizer == "optim") {
    result = stats::optim(eta, objective, gradient, method = control$method, control = options)
  } else {
    arguments = list(logLik = function(value) evaluator$evaluate(value)$loglik,
      grad = function(value) evaluator$evaluate(value)$optimizer_score,
      start = eta, method = control$method, control = options)
    # maxLik's optim-based methods take parscale as an argument, not MaxControl.
    if ("parscale" %in% names(options)) {
      if (toupper(control$method) %in% c("NR", "NEWTON-RAPHSON", "BHHH", "BFGSR", "BFGS-R"))
        stop("parscale is supported by optim-based methods, not native maxLik methods", call. = FALSE)
      arguments$parscale = options$parscale
      arguments$control$parscale = NULL
    }
    # OPG requires no final Hessian. Explicit observed inference is computed
    # below from the analytic score; avoid BFGS's duplicate final calculation.
    # A Newton optimizer can still require curvature during its iterations.
    if (control$vcov == 'opg' || toupper(control$method) == 'BFGS')
      arguments$finalHessian = FALSE
    raw = do.call(maxlik_namespace()$maxLik, arguments)
    result = list(par = raw$estimate, convergence = raw$code,
                   counts = c("function" = unname(raw$iterations), gradient = NA_integer_),
                   message = raw$message, raw = raw)
  }
  result
}

# maxLik is suggested rather than imported, so it is loaded only when a fit
# asks for it; `installed` lets tests exercise the message.
maxlik_namespace = function(installed = requireNamespace("maxLik", quietly = TRUE)) {
  if (!installed)
    stop("optimizer = \"maxLik\" needs the maxLik package; install it or use optimizer = \"optim\"",
         call. = FALSE)
  asNamespace("maxLik")
}

# maxLik's Newton and optim-based methods assign different meanings to codes.
optimizer_status = function(code, control, options, n_parameters) {
  if (!n_parameters) return("all_fixed")
  limit = if (control$optimizer == "optim") options$maxit else options$iterlim
  if (identical(as.integer(limit), 0L)) return("not_run")
  # Both backends return code 0 for SANN whatever its state; it has no stopping test.
  if (toupper(control$method) == "SANN") return("no_convergence_test")
  newton = control$optimizer == "maxLik" && toupper(control$method) %in%
    c("NR", "NEWTON-RAPHSON", "BHHH", "BFGSR", "BFGS-R")
  if (code %in% if (newton) c(1L, 2L, 8L) else 0L) "converged" else "not_converged"
}

# 3. Fixed optimizer scaling and stopping controls ---------------------------

# Unit reference for each coordinate's panel-score column norm, NA except for
# regression coefficients: the design scale times sqrt(panels), divided by the
# outcome scale (root mean square of the selected outcomes) for Linear outcome
# coefficients, whose scores carry 1/y units. Ratios of score norms to it are
# unchanged by rescaled regressors or Linear outcomes.
score_references = function(prepared, names) {
  scales = if (!is.null(prepared$x)) design_scales(prepared)
  if (!length(scales)) return(rep(NA_real_, length(names)))
  reference = unname(scales[names])*sqrt(prepared$n_panels)
  y = if (prepared$family == "gaussian") prepared$y[prepared$selected]
  y = y[is.finite(y) & y != 0]
  outcome = startsWith(names, "outcome.")
  if (length(y)) reference[outcome] = reference[outcome]/sqrt(mean(y^2))
  reference
}

optimizer_options = function(control, prepared = NULL, info = NULL, score = NULL) {
  options = list()
  if (control$optimizer == "optim") {
    options$maxit = control$maxit
    options$trace = control$trace
    if (control$method == "L-BFGS-B") options$factr = control$factr else
      options$reltol = control$reltol
  } else {
    options$iterlim = control$maxit
    options$reltol = control$reltol
    options$printLevel = control$trace
  }
  if (!is.null(control$parscale)) options$parscale = control$parscale
  # Joint fixed binary fits can propose costly extreme-tail trials with BFGS's
  # identity initial metric. Scale by the number of independent panels, paired
  # with tighter stopping. Preserve explicit controls and other model routes.
  if (!is.null(prepared) && prepared$model_id == 'probit_re_probit_re' &&
      control$quadrature == 'gh' &&
      control$optimizer == 'optim' && control$method == 'BFGS' &&
      'cor_error' %in% info$name[!info$fixed] &&
      isFALSE(attr(control, 'reltol_supplied', exact = TRUE)) &&
      identical(control$reltol, sqrt(.Machine$double.eps)) &&
      is.null(control$parscale)) {
    options$parscale = rep(1/sqrt(prepared$n_panels), sum(!info$fixed))
    options$reltol = 1e-12
  }
  if (!"parscale" %in% names(options) && toupper(control$method) == "BFGS" &&
      !is.null(score) && ncol(score)) {
    # A fixed diagonal metric from the already-cached initial panel scores.
    # Each scale is the reciprocal root sum of squared unit scores, the
    # reciprocal square root of the OPG diagonal. Unit scores can retain
    # variation when their total is zero at an optimal start. Coefficient
    # scales follow their units, so rescaling a regressor or a Linear outcome
    # rescales the BFGS path. A structurally zero coefficient score (as at a
    # symmetric zero start), below 1e-6 of its unit reference, and the
    # unit-free log-SD and atanh-correlation coordinates keep the cap at one
    # against unbounded weak-direction steps. The lower bound avoids zero
    # scales on overflow.
    norm = sqrt(colSums(score^2))
    reference = score_references(prepared, colnames(score) %||% character(length(norm)))
    coefficient = !is.na(reference) & norm > 1e-6*reference
    options$parscale = pmax(sqrt(.Machine$double.eps), 1/unname(ifelse(coefficient, norm, pmax(1, norm))))
  }
  if (!is.null(info) && "parscale" %in% names(options)) {
    scale = options$parscale; free = info$name[!info$fixed]
    if (!is.numeric(scale) || is.complex(scale) || !length(scale) ||
        any(!is.finite(scale) | scale <= 0) || !length(scale) %in% c(1L, length(free)))
      stop("parscale must have one positive finite value or one per free parameter", call. = FALSE)
    if (length(scale) > 1L && !is.null(names(scale))) {
      if (anyDuplicated(names(scale)) || !setequal(names(scale), free))
        stop("Named parscale values must match the free parameter names", call. = FALSE)
      scale = scale[free]
    }
    options$parscale = rep_len(unname(scale), length(free))
  }
  options
}
