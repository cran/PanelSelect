# Exact Gaussian random-intercept auxiliary for selected outcome rows.
# This supplies starts only; it is not a public standalone Gaussian model.

# 1. Preserve the prepared outcome design and natural parameter restrictions
linear_re_initial_fit = function(prepared, start, fixed = numeric(), control = panel_control()) {
  parameters = c(if (ncol(prepared$x)) paste0("outcome.", colnames(prepared$x)),
    "sd_re_outcome", "sd_error")
  info = data.frame(name = parameters,
    type = ifelse(startsWith(parameters, "sd_"), "sd", "coefficient"),
    fixed = parameters %in% names(fixed), stringsAsFactors = FALSE)
  start = validate_named_parameters(start, info, "Gaussian auxiliary start", allow_zero = TRUE)
  fixed = validate_named_parameters(fixed, info, "Gaussian auxiliary fixed", allow_zero = TRUE)
  if (!setequal(names(start), parameters))
    stop("Gaussian auxiliary starts must contain the complete outcome parameter block", call. = FALSE)
  if (any(start[names(fixed)] != fixed))
    stop("A Gaussian auxiliary start conflicts with a fixed parameter", call. = FALSE)
  theta = start[parameters]
  free = !info$fixed
  if (theta[["sd_error"]] <= 0 || any(theta[info$type == "sd" & free] <= 0))
    stop("Gaussian auxiliary residual and free standard deviations must be positive", call. = FALSE)
  if (!any(prepared$selected))
    stop("Gaussian outcome initialization needs selected rows", call. = FALSE)

  # Selection is absent from this auxiliary likelihood. The shared marginal
  # uses only selected x, y and offsets; panels without outcomes contribute zero.
  # Keep its own preparation, without mutating any shared integration memo.
  prepared$integration_cache = NULL
  prepared = prepare_model_data(prepared)

  # 2. Exact marginal likelihood and scores in free optimizer coordinates
  # For m selected rows, covariance is sigma^2 I + lambda^2 11'. Reuse the
  # centred Gaussian calculation and its natural scores from the model file.
  cached_eta = cached_value = NULL
  evaluate = function(eta, trial = TRUE) {
    if (!identical(eta, cached_eta)) {
      at = parameter_from_optimizer(eta, theta, info, check = FALSE)
      valid = all(is.finite(at)) && all(at[info$type == "sd" & free] > 0)
      value = if (valid) linear_outcome_marginal(at, prepared, conditional = FALSE)
      valid = valid && all(is.finite(value$value)) && all(is.finite(value$gradient[, free, drop = FALSE]))
      if (valid) {
        score = value$gradient[, free, drop = FALSE]
        colnames(score) = parameters[free]
        score = sweep(score, 2L, parameter_jacobian(at, info, check = FALSE)[free], `*`)
        valid = all(is.finite(score))
      }
      cached_value <<- list(theta = at, valid = valid,
        loglik = if (valid) value$value else rep(-Inf, prepared$n_panels),
        optimizer_score = if (valid) score else matrix(0, prepared$n_panels, sum(free)))
      cached_eta <<- eta
    }
    # Reject unrepresentable line-search trials, but never accept an invalid
    # initial or final evaluation as an auxiliary fit.
    if (!trial && !cached_value$valid)
      stop("Gaussian auxiliary likelihood or score is not finite", call. = FALSE)
    cached_value
  }

  # 3. Optimize once; the caller decides whether to retain these starts
  # No numerical integration, observed information or covariance is requested.
  eta = parameter_to_optimizer(theta, info)
  initial = evaluate(eta, trial = FALSE)
  options = optimizer_options(control, prepared, info, initial$optimizer_score)
  optimizer = optimize_model(eta, list(evaluate = evaluate), control, options)
  final = evaluate(optimizer$par, trial = FALSE)
  optimizer$status = optimizer_status(optimizer$convergence, control, options, length(eta))
  optimizer$converged = optimizer$status %in% c("converged", "all_fixed")
  list(coefficients = final$theta, optimizer = optimizer, loglik = sum(final$loglik))
}
