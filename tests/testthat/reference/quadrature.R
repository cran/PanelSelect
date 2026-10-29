# Frozen 3ffd4c7 forward adaptive rule and cold mode/geometry calculations.
# These are regression oracles, not independent references: they reuse the
# production model kernels (make_panel_kernel, make_model_context) and freeze an
# earlier assembly of them. Agreement shows the current engine reproduces that
# assembly; independent values come from reference/math.R and the fixtures.

reference_integrate_prepared = function(theta, prepared, control, fixed, score = TRUE,
                              posterior = FALSE) {
  control = prepare_rule_control(control)
  q = model_dimension(theta, prepared, control, fixed)
  n = prepared$n_panels; p = length(theta)
  differentiated = score || posterior
  adaptive = q > 0 && (control$quadrature == 'agh' ||
                       !is.null(control$fixed_geometry))
  context = make_model_context(theta, prepared, control, fixed,
    indices_only = !prepared$has_selection && !adaptive && !posterior)
  geometry = NULL
  if (adaptive) {
    geometry = control$fixed_geometry
    if (is.null(geometry)) {
      geometry = panel_geometry(theta, prepared, control, fixed, differentiated, context)
    } else {
      if (!identical(dim(geometry$mode), c(n, q)))
        integration_error('Frozen integration geometry does not match the panels')
      geometry$mode_gradient = array(0, c(n, q, p))
      geometry$transform_gradient = array(0, c(n, q, q, p))
      geometry$jacobian_gradient = matrix(0, n, p)
    }
  }
  # The exact Gaussian outcome marginal is common to every outer node. Mode
  # search above still sees the full kernel, including its acceptance tolerance.
  marginal = if (!posterior) context$marginal else NULL
  if (!is.null(marginal)) context$marginal = NULL
  # With both Gaussian correlations fixed to zero, the conditional selection
  # index is standalone probit. Geometry above retains the full marginal during
  # mode acceptance; only the outer value/first-score calculation is simplified.
  kernel_data = prepared
  if (!posterior && context$gaussian && context$tau_zero &&
      is_fixed_zero('cor_re', theta, fixed)) {
    kernel_data = prepared$integration_cache$selection %||%
      prepare_model_data(selection_component(prepared))
    context = make_model_context(theta, kernel_data, control, fixed)
  }
  direct = !kernel_data$has_selection && !adaptive && !posterior
  kernel = if (!direct) make_panel_kernel(theta, kernel_data, control, fixed,
    order = as.integer(differentiated), parameters = differentiated,
    latent = adaptive, context = context)
  # The rule cache keeps one-dimensional rules per effect (selection, outcome).
  rule = if (q == 1L) control$.rules$outer[[2L]][[single_effect(theta, prepared)]] else
    control$.rules$outer[[q + 1L]]
  count = nrow(rule$nodes)
  batch = quadrature_block_size(prepared, q, p, control$observation_nodes, differentiated, count)
  state = NULL
  saved = if (posterior) lapply(seq_len(n), function(i) list(
    nodes = matrix(0, count, 2L, dimnames = list(NULL, c('u', 'v'))),
    node_gradient = array(0, c(count, 2L, p), dimnames = list(NULL, c('u', 'v'), names(theta))),
    log_weights = numeric(count),
    weight_score = matrix(0, count, p, dimnames = list(NULL, names(theta))))) else NULL

  for (first in seq.int(1L, count, by = batch)) {
    indices = seq.int(first, min(count, first + batch - 1L))
    base = rule$nodes[indices, , drop = FALSE]
    size = length(indices)
    if (direct) {
      reduced = integrate_index_block(theta, kernel_data, context, base,
        rule$log_weights[indices], differentiated)
      state = merge_quadrature(state, reduced$loglik, reduced$score)
      next
    }
    if (adaptive) {
      mapped = affine_batch(base, geometry$mode, geometry$transform,
        if (differentiated) geometry$mode_gradient else NULL,
        if (differentiated) geometry$transform_gradient else NULL)
      z = mapped$nodes
    } else z = base[rep(seq_len(size), each = n), , drop = FALSE]
    result = kernel$evaluate(z, prior = adaptive)
    value = result$value + rep(rule$log_weights[indices], each = n)
    derivative = if (differentiated)
      result$gradient[, kernel$parameter_columns, drop = FALSE] else matrix(0, n*size, p)
    if (adaptive) {
      value = value + rep(geometry$log_jacobian, size) +
        rep(rowSums(base^2)/2 + q/2*log(2*pi), each = n)
      if (differentiated) derivative = derivative +
        geometry$jacobian_gradient[rep(seq_len(n), size), , drop = FALSE] +
        row_contract(result$gradient[, seq_len(q), drop = FALSE], mapped$gradient)
    }
    reduced = logsum_score_cpp(value, derivative, n)
    state = merge_quadrature(state, reduced[, 1L], reduced[, -1L, drop = FALSE])

    if (posterior) {
      effects = kernel$effects(z)
      for (a in 1:2) {
        effect = effects[[a]]
        effect_gradient = effect$gradient[, kernel$parameter_columns, drop = FALSE]
        if (adaptive) effect_gradient = effect_gradient +
          row_contract(effect$gradient[, seq_len(q), drop = FALSE], mapped$gradient)
        for (i in seq_len(n)) {
          rows = i + n*(seq_len(size) - 1L)
          saved[[i]]$nodes[indices, a] = effect$value[rows]
          saved[[i]]$node_gradient[indices, a, ] = effect_gradient[rows, , drop = FALSE]
        }
      }
      for (i in seq_len(n)) {
        rows = i + n*(seq_len(size) - 1L)
        saved[[i]]$log_weights[indices] = value[rows]
        saved[[i]]$weight_score[indices, ] = derivative[rows, , drop = FALSE]
      }
    }
  }
  if (!is.null(marginal)) {
    state$loglik = state$loglik + marginal$value
    if (differentiated) state$score = state$score + marginal$gradient
  }
  if (any(!is.finite(state$loglik)) || (score && any(!is.finite(state$score))))
    integration_error('Likelihood or score is not finite at this parameter vector')
  colnames(state$score) = names(theta)
  state$diagnostics = list(method = if (adaptive) 'agh' else 'gh', dimension = q,
    nodes = control$nodes, observation_nodes = control$observation_nodes,
    mode_residual = if (adaptive) geometry$mode_residual else NULL,
    order_checked = FALSE)
  if (posterior) {
    for (i in seq_len(n)) saved[[i]]$log_weights = saved[[i]]$log_weights - state$loglik[i]
    names(saved) = as.character(prepared$id_values)
    state$posterior = saved
  }
  state
}


reference_panel_geometry = function(theta, prepared, control, fixed, derivatives = TRUE,
                          context = NULL) {
  if (is.null(context)) context = make_model_context(theta, prepared, control, fixed)
  kernel = make_panel_kernel(theta, prepared, control, fixed,
    order = 2L, parameters = FALSE, latent = TRUE, context = context)
  values = make_panel_kernel(theta, prepared, control, fixed,
    order = 0L, parameters = FALSE, latent = TRUE, context = context)
  q = kernel$q; n = prepared$n_panels; p = length(theta)
  mode = matrix(0, n, q)
  current = kernel$evaluate(mode, prior = TRUE)
  done = rep(FALSE, n)
  iterations = 0L
  for (iteration in seq_len(control$mode_maxit)) {
    iterations = iteration
    if (any(!is.finite(current$value)) || any(!is.finite(current$gradient)) ||
        any(!is.finite(current$hessian)))
      integration_error("Adaptive mode has nonfinite likelihood derivatives")
    gradient = current$gradient
    done = row_max_abs(gradient) <= control$mode_tol*(1 + row_max_abs(mode))
    if (all(done)) break
    direction = matrix(0, n, q)
    for (i in which(!done)) {
      information = -matrix(current$hessian[i, , ], q, q)
      factor = tryCatch(chol(information), error = function(e) NULL)
      step = if (is.null(factor)) gradient[i, ] else
        drop(backsolve(factor, forwardsolve(t(factor), gradient[i, ])))
      direction[i, ] = step/max(1, sqrt(sum(step^2))/5)
    }
    if (any(!is.finite(direction))) integration_error("Adaptive mode step is not finite")
    multiplier = rep(1, n)
    accepted = done
    for (backtrack in seq_len(40L)) {
      proposal = mode + direction*multiplier
      next_value = values$evaluate(proposal, prior = TRUE)
      improvement = next_value$value - current$value
      required = 1e-4*multiplier*rowSums(gradient*direction)
      admissible = is.finite(next_value$value) & is.finite(improvement) &
        is.finite(required) & improvement >= required - 1e-13*(1 + abs(current$value))
      admissible[is.na(admissible)] = FALSE
      accepted = accepted | admissible
      if (all(accepted)) break
      multiplier[!accepted] = multiplier[!accepted]/2
    }
    if (!all(accepted)) integration_error("Adaptive mode line search failed")
    mode = proposal
    current = kernel$evaluate(mode, prior = TRUE)
  }
  residual = row_max_abs(current$gradient)
  if (any(!is.finite(residual)) ||
      any(residual > 10*control$mode_tol*(1 + row_max_abs(mode))))
    integration_error("Adaptive mode did not satisfy its gradient tolerance")

  full = if (derivatives) make_panel_kernel(theta, prepared, control, fixed,
    order = 3L, parameters = TRUE, latent = TRUE, context = context)$evaluate(mode, prior = TRUE) else current
  transform = array(0, c(n, q, q))
  mode_gradient = array(0, c(n, q, p))
  transform_gradient = array(0, c(n, q, q, p))
  log_jacobian = numeric(n)
  jacobian_gradient = matrix(0, n, p)
  identity = diag(q)
  upper = upper.tri(identity)
  diagonal = seq.int(1L, q*q, by = q+1L)
  for (i in seq_len(n)) {
    information = -matrix(full$hessian[i, , seq_len(q)], q, q)
    factor = tryCatch(chol(information), error = function(e) NULL)
    if (is.null(factor)) integration_error("Adaptive mode curvature is not positive definite")
    lower = t(factor)
    inverse = chol2inv(factor)
    C = backsolve(factor, identity)
    transform[i, , ] = C
    log_jacobian[i] = -sum(log(diag(factor)))
    if (!derivatives) next
    inverse_lower_transpose = t(forwardsolve(lower, identity))
    dm = inverse %*% matrix(full$hessian[i, , q + seq_len(p)], q, p)
    mode_gradient[i, , ] = dm
    for (j in seq_len(p)) {
      dH = -matrix(full$third[i, , , q + j], q, q)
      for (a in seq_len(q))
        dH = dH - matrix(full$third[i, , , a], q, q)*dm[a, j]
      K = forwardsolve(lower, dH) %*% inverse_lower_transpose
      K[upper] = 0
      K[diagonal] = K[diagonal]/2
      dL = lower %*% K
      transform_gradient[i, , , j] = -C %*% t(dL) %*% C
      jacobian_gradient[i, j] = -.5*sum(inverse*t(dH))
    }
  }
  list(mode = mode, transform = transform, mode_gradient = mode_gradient,
       transform_gradient = transform_gradient, log_jacobian = log_jacobian,
       jacobian_gradient = jacobian_gradient, mode_residual = residual,
       iterations = iterations)
}

