# Integration over person effects.
#
# Start at integrate_prepared: it calls the model card's likelihood entry,
# then checks and describes the per-unit result. Most entries use
# integrate_person_effects, which reduces bounded batches of node values and
# scores from make_panel_kernel. Model entries may instead supply collapsed
# fixed-rule blocks (integrate_index_block, integrate_selection_fixed).
#
# Fixed GH differentiates a fixed finite sum. Adaptive GH also differentiates
# its mode, scale and log-volume correction: panel_geometry/panel_mode build
# that geometry and adaptive_score assembles the moving-node terms. The local
# likelihoods and natural first-score chains remain visible in model files;
# model_assembly.R supplies shared sums and higher mixed derivatives.
# posterior_rule retains nodes and their derivatives for history predictions.

# Byte bounds for temporary derivative and quadrature buffers (16 MiB) and for
# staged outer GL batches and workspaces (64 MiB); batches stay below them.
buffer_bytes = 16*1024^2
staged_buffer_bytes = 64*1024^2
# Refinement passes of one staged mesh before it is reported unresolved.
staged_refinement_passes = 128L

normal_rule = function(nodes, q) {
  if (!q) return(list(nodes = matrix(numeric(), 1L, 0L), log_weights = 0))
  nodes = rep_len(nodes, q)
  rules = lapply(nodes, function(n) statmod::gauss.quad(n, "hermite"))
  index = as.matrix(expand.grid(lapply(nodes, seq_len)))
  points = matrix(0, nrow(index), q)
  weights = rep(0, nrow(index))
  for (j in seq_len(q)) {
    points[, j] = sqrt(2)*rules[[j]]$nodes[index[, j]]
    weights = weights + log(rules[[j]]$weights[index[, j]]) - .5*log(pi)
  }
  list(nodes = points, log_weights = weights)
}

# `dimensions` is the largest latent dimension the caller can evaluate: two for
# panel models, zero for cross-sectional ones, whose rules are never used.
prepare_rule_control = function(control, dimensions = 2L) {
  rules = control$.rules
  observation_quadrature = control$observation_quadrature %||% 'gh'
  # A saved control with a removed setting stops with its refit, and rule names
  # are checked, so no name from another vocabulary falls through to fixed GH.
  stop_removed_settings(control)
  if (!isTRUE(control$quadrature %in% c('agh', 'gh', 'gl')) ||
      !isTRUE(observation_quadrature %in% c('agh', 'gh', 'gl')))
    stop('Internal error: unknown integration rule ', control$quadrature, ' or ',
         observation_quadrature, call. = FALSE)
  outer_gh = !uses_gl_quadrature(control)
  if (is.null(rules) || !identical(rules$layout, 2L) || !identical(rules$nodes, control$nodes) ||
      !identical(rules$observation_nodes, control$observation_nodes) ||
      !identical(rules$observation_quadrature %||% 'gh', observation_quadrature) ||
      !identical(rules$outer_gh %||% TRUE, outer_gh) || (rules$dimensions %||% 2L) < dimensions) {
    # Outer rules by dimension; one dimension has a rule for each effect.
    orders = rep_len(control$nodes, 2L)
    control$.rules = list(layout = 2L, nodes = control$nodes,
      observation_nodes = control$observation_nodes,
      observation_quadrature = observation_quadrature, outer_gh = outer_gh, dimensions = dimensions,
      outer = if (outer_gh) c(list(normal_rule(orders, 0L)),
        if (dimensions >= 1L) list(list(selection = normal_rule(orders[1L], 1L),
          outcome = normal_rule(orders[2L], 1L))),
        if (dimensions >= 2L) list(normal_rule(orders, 2L))) else NULL,
      inner = if (observation_quadrature %in% c('gh', 'agh'))
        normal_rule(control$observation_nodes, 1L) else NULL)
  }
  control
}

quadrature_block_size = function(prepared, q, p, observation_nodes, score, count,
                                contracted = FALSE) {
  local = if (prepared$family == 'pln') 4L else 3L
  inner = if (prepared$family == 'pln') observation_nodes else 1L
  row_entries = prepared$n_rows*(8 + 4*local*inner + if (score) 2*p else 0)
  geometry_entries = if (contracted) q + q*(q+1L)/2 else 2*q*p
  panel_entries = prepared$n_panels*(4 + 4*p + q + if (score) geometry_entries else 0)
  max(1L, min(count, floor(buffer_bytes/(8*max(1, row_entries + panel_entries)))))
}

integration_error = function(message) {
  stop(structure(list(message = message, call = NULL),
                 class = c("panelselect_integration_error", "error", "condition")))
}

merge_quadrature = function(state, value, score) {
  if (is.null(state)) return(list(loglik = value, score = score))
  high = pmax(state$loglik, value)
  left = exp(state$loglik - high); right = exp(value - high)
  absent = !is.finite(high) & high < 0
  left[absent] = right[absent] = 0
  total = left + right
  a = ifelse(total > 0, left/total, 0)
  b = ifelse(total > 0, right/total, 0)
  old_score = state$score; new_score = score
  old_score[a == 0, ] = 0; new_score[b == 0, ] = 0
  list(loglik = ifelse(absent, -Inf, high + log(total)),
       score = old_score*a + new_score*b)
}

# Fixed GH for a standalone model: B_tk = b_t + lambda*z_k at every node, the
# model's outcome rows, then node partials averaged under the posterior weights
# before the design is expanded (standalone_fixed_score).
integrate_index_block = function(theta, prepared, context, base, log_weights, score) {
  block = context$data$active; count = length(block$rows)
  groups = prepared$n_panels; nodes = nrow(base)
  gradient = matrix(0, groups, length(theta))
  value = matrix(0, groups, nodes)
  if (count) {
    index = rep(context$b[block$rows], nodes)
    if (context$q) index = index + rep(context$lambda*base[, 1L], each = count)
    local = context$definition$outcome_rows(index, rep(block$y, nodes),
      rep(context$data$log_factorial[block$rows], nodes), context, as.integer(score))
    value = group_sum(matrix(local$value, count, nodes), block$ptr)
  }
  mixture = logsum_weights(as.vector(value) + rep(log_weights, each = groups), groups)
  if (count && score) {
    # array() keeps the leading count*nodes values: column 1 (the B partial)
    # of the column-major gradient, in one copy.
    first = array(local$gradient, c(count, nodes))
    product = group_weighted_design(block$x, first, block$ptr, mixture$weights)
    sums = list(B = if (context$q) group_sum(first, block$ptr),
      gamma = if ('error' %in% local$coordinates)
        group_sum(matrix(local$gradient[, 2L], count, nodes), block$ptr))
    gradient = standalone_fixed_score(gradient,
      product[, seq_along(context$ix), drop = FALSE], sums,
      if (context$q) base[, 1L] else NULL, mixture$weights, context$score_columns)
  }
  list(loglik = mixture$value, score = gradient)
}

# Posterior-weighted natural scores of a standalone model under a fixed rule:
# beta takes the mixture-weighted design products; lambda E[v sum_t l_B]; and
# gamma (PLN) E[sum_t l_gamma]. sums and weights are panel x node. Exactly zero
# weights mask their nodes, so an infinite partial at an underflowed node
# cannot turn the score into NaN.
standalone_fixed_score = function(score, products, sums, nodes, weights, columns) {
  if (length(columns$beta)) score[, columns$beta] = products
  if (length(columns$lambda) && !is.null(nodes)) {
    d_lambda = sweep(sums$B, 2L, nodes, '*')
    d_lambda[weights == 0] = 0
    score[, columns$lambda] = rowSums(d_lambda*weights)
  }
  if (length(columns$gamma) && !is.null(sums$gamma)) {
    d_gamma = sums$gamma
    d_gamma[weights == 0] = 0
    score[, columns$gamma] = rowSums(d_gamma*weights)
  }
  score
}

# With independent observation errors, the mass is Phi((2d-1)A)*f(y|B)^d.
# Condition on u, integrate the outcome nodes, then integrate u. This is the
# same two-dimensional GH rule, including when free rho happens to equal zero.
integrate_selection_fixed = function(theta, prepared, context, rule, score,
                                     batch = NULL) {
  n = prepared$n_panels; p = length(theta)
  selection = context$data$active; outcome = context$data$selected
  count = length(outcome$rows)
  # The tensor rule stores u fastest. Recover its two marginal log weights
  # from the cached rule, without regenerating GH nodes at every evaluation.
  u_nodes = unique(rule$nodes[, 1L]); z = unique(rule$nodes[, 2L])
  staged_record_outer_measure(context, prepared, rule$nodes, rule$log_weights, shared = TRUE)
  log_weights = matrix(rule$log_weights, length(u_nodes), length(z))
  # Extreme GH weights may underflow. Recover factors from a positive-mass
  # row/column, so zero-mass edge nodes keep log weight -Inf rather than NaN.
  pivot = which(is.finite(log_weights), arr.ind = TRUE)[1L, ]
  log_u = log_weights[, pivot[2L]]
  log_u = log_u-logsum_weights(log_u, 1L)$value
  log_z = log_weights[pivot[1L], ]-log_u[pivot[1L]]
  # Bound additional observation/node buffers to about 16 MiB. The required
  # panel x parameter result is independent of the number of quadrature nodes.
  if (is.null(batch)) batch = max(1L, min(length(z),
    floor(buffer_bytes/(8*max(1, 8*count+4*n)))))
  rho = theta[['cor_re']]; c_rho = sqrt(1-rho^2)
  state = NULL
  for (k in seq_along(u_nodes)) {
    u = u_nodes[k]
    # log Phi(s_t A_t), s_t = 2d_t - 1: the selection factor at this u node.
    A = context$a[selection$rows]+context$delta*u
    local = local_univariate_probit(A, 2*prepared$d[selection$rows]-1, 'a', as.integer(score))
    log_selection = drop(group_sum(matrix(local$value), selection$ptr))
    # The single-column gradient is already the rows x 1 matrix of A partials.
    selection_score = if (score) group_weighted_design(selection$w,
      local$gradient, selection$ptr) else NULL
    conditional = NULL
    for (first in seq.int(1L, length(z), by = batch)) {
      indices = seq.int(first, min(length(z), first+batch-1L))
      v = rho*u+c_rho*z[indices]; size = length(indices)
      B = rep(context$b[outcome$rows], size)+rep(context$lambda*v, each = count)
      local = context$definition$outcome_rows(B, rep(outcome$y, size),
        rep(context$data$log_factorial[outcome$rows], size), context, as.integer(score))
      log_outcome = group_sum(matrix(local$value, count, size), outcome$ptr)
      mixture = logsum_weights(as.vector(log_outcome)+rep(log_z[indices], each = n), n)
      gradient = matrix(0, n, p)
      if (score) {
        # Column 1 (the B partial) as rows x nodes; see integrate_index_block.
        d_B = array(local$gradient, c(count, size))
        products = group_weighted_design(outcome$x, d_B, outcome$ptr, mixture$weights)
        gradient = selection_fixed_score(gradient,
          products[, seq_along(context$ix), drop = FALSE],
          list(B = group_sum(d_B, outcome$ptr),
            gamma = if ('error' %in% local$coordinates)
              group_sum(matrix(local$gradient[, 2L], count, size), outcome$ptr)),
          list(v = v, v_rho = u-rho/c_rho*z[indices], lambda = context$lambda),
          mixture$weights, context$score_columns)
      }
      conditional = merge_quadrature(conditional, mixture$value, gradient)
    }
    # Selection scores at this u node: g_alpha = sum w l_A, g_delta = u sum l_A.
    if (score) conditional$score = chain_score(conditional$score,
      list(list(alpha = 'w', delta = u)), list(selection_score[, length(context$iw)+1L]),
      list(selection_score), list(columns = context$score_columns), NULL)
    state = merge_quadrature(state, log_selection+conditional$loglik+log_u[k],
      conditional$score)
  }
  state
}

# Conditional outcome scores after integrating z at a fixed selection node u,
# for the outcome index B = b + lambda*v with v = rho*u + sqrt(1-rho^2)*z:
# g_beta, g_lambda = E[v sum_t l_B], g_rho = lambda E[(dv/drho) sum_t l_B] and
# g_gamma (PLN). sums and weights are panel x outcome-node; products has one
# column per beta. integrate_selection_fixed adds the selection alpha/delta blocks.
selection_fixed_score = function(score, products, sums, nodes, weights, columns) {
  if (length(columns$beta)) score[, columns$beta] = products
  d_B = sums$B
  d_B[weights == 0] = 0
  weighted = d_B*weights
  if (length(columns$lambda)) score[, columns$lambda] = drop(weighted %*% nodes$v)
  if (length(columns$rho)) score[, columns$rho] =
    nodes$lambda*drop(weighted %*% nodes$v_rho)
  if (length(columns$gamma) && !is.null(sums$gamma)) {
    d_gamma = sums$gamma
    d_gamma[weights == 0] = 0
    score[, columns$gamma] = rowSums(d_gamma*weights)
  }
  # Preserve E[v*d_B] and lambda*E[v_rho*d_B], not products of expectations.
  score
}

# One parameter point: the model's own likelihood entry (card$likelihood,
# integrate_person_effects by default) returns the log likelihood and natural
# score of every independent unit; they are checked and described here.
integrate_prepared = function(theta, prepared, control, fixed, score = TRUE,
                              posterior = FALSE) {
  q = model_dimension(theta, prepared, control, fixed)
  control = prepare_rule_control(control, q)
  model = resolve_model(prepared$family, prepared$has_selection, prepared$structure,
    prepared$model_id)
  state = model$likelihood(theta, prepared, control, fixed, q, score, posterior)
  result = list(loglik = state$loglik, score = state$score)
  if (any(!is.finite(result$loglik)) || (score && any(!is.finite(result$score))))
    integration_error('Likelihood or score is not finite at this parameter vector')
  colnames(result$score) = names(theta)
  # A likelihood with no remaining latent effect uses no quadrature rule.
  adaptive = q > 0 && control$quadrature == 'agh'
  remaining = if (q == 1L) single_effect(theta, prepared)
  result$diagnostics = list(method = if (q == 0L) 'none' else if (adaptive) 'agh' else 'gh',
    dimension = q,
    nodes = if (q == 1L) rep_len(control$nodes, 2L)[[if (remaining == 'selection') 1L else 2L]] else
      control$nodes, observation_nodes = control$observation_nodes,
    mode_residual = if (adaptive) state$mode_residual else NULL,
    order_checked = FALSE)
  if (posterior) {
    saved = state$posterior
    for (i in seq_len(prepared$n_panels))
      saved[[i]]$log_weights = saved[[i]]$log_weights - result$loglik[i]
    names(saved) = as.character(prepared$id_values)
    result$posterior = saved
  }
  result
}

# Integrates every unit's panel kernel over its active person effects:
# log L_i = log sum_k w_k exp(r_i(z_k)) at fixed GH nodes z_k or, for adaptive
# GH, at z = mode + C x with the log-volume correction. Fixed-rule scores are
# posterior-weighted node scores; adaptive scores also include derivatives of
# the moving nodes and log volume (adaptive_score).
# m is the model point (make_model_context). block replaces the kernel by a
# model's collapsed fixed-rule block (integrate_index_block,
# probit_re_linear_re_fixed_block); kernel_point and kernel_data integrate
# another factor than m (Linear selection with both correlations fixed at 0).
# Returns loglik and score per unit, the adaptive mode residual and, for history
# predictions, the posterior nodes before normalization.
integrate_person_effects = function(theta, prepared, control, fixed, q, score,
    posterior = FALSE, m = make_model_context(theta, prepared, control, fixed),
    block = NULL, kernel_point = NULL, kernel_data = prepared) {
  n = prepared$n_panels; p = length(theta)
  differentiated = score || posterior
  adaptive = q > 0 && control$quadrature == 'agh'
  contracted = adaptive && differentiated && !posterior
  geometry = if (adaptive) panel_geometry(theta, prepared, control, fixed, differentiated, m,
    contracted = contracted)
  # The exact Gaussian outcome marginal is common to every outer node. Mode
  # search above still sees the full kernel, including its acceptance tolerance.
  marginal = if (!posterior) m$marginal else NULL
  if (!is.null(marginal)) m$marginal = NULL
  if (is.null(kernel_point)) kernel_point = m
  kernel = if (is.null(block)) make_panel_kernel(theta, kernel_data, control, fixed,
    order = as.integer(differentiated), parameters = differentiated,
    latent = adaptive, context = kernel_point)
  remaining = if (q == 1L) single_effect(theta, prepared)
  rule = if (q == 1L) control$.rules$outer[[2L]][[remaining]] else control$.rules$outer[[q + 1L]]
  count = nrow(rule$nodes)
  batch = quadrature_block_size(prepared, q, p, control$observation_nodes, differentiated, count,
    contracted = contracted)
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
    if (!is.null(block)) {
      # A collapsed fixed rule integrates this batch of nodes itself.
      staged_record_outer_measure(kernel_point, prepared, base, rule$log_weights[indices], shared = TRUE)
      reduced = block(theta, kernel_data, kernel_point, base, rule$log_weights[indices],
        differentiated)
      state = merge_quadrature(state, reduced$loglik, reduced$score)
      next
    }
    if (adaptive) {
      mapped = affine_batch(base, geometry$mode, geometry$transform,
        if (differentiated && !contracted) geometry$mode_gradient else NULL,
        if (differentiated && !contracted) geometry$transform_gradient else NULL)
      z = mapped$nodes
    } else z = base[rep(seq_len(size), each = n), , drop = FALSE]
    measure_weights = rep(rule$log_weights[indices], each = n)
    if (adaptive) measure_weights = measure_weights+rep(geometry$log_jacobian, size)+
      rep(rowSums(base^2)/2, each = n)-rowSums(z^2)/2
    staged_record_outer_measure(kernel_point, prepared, z, measure_weights)
    result = kernel$evaluate(z, prior = adaptive)
    value = result$value + rep(rule$log_weights[indices], each = n)
    derivative = if (differentiated)
      result$gradient[, kernel$parameter_columns, drop = FALSE] else matrix(0, n*size, p)
    if (adaptive) {
      value = value + rep(geometry$log_jacobian, size) +
        rep(rowSums(base^2)/2 + q/2*log(2*pi), each = n)
      if (contracted) {
        latent = result$gradient[, seq_len(q), drop = FALSE]
        moments = latent[, 1L]*rep(base[, 1L], each = n)
        if (q == 2L) moments = cbind(moments,
          latent[, 1L]*rep(base[, 2L], each = n), latent[, 2L]*rep(base[, 2L], each = n))
        derivative = cbind(derivative, latent, moments)
      } else if (differentiated) derivative = derivative +
          geometry$jacobian_gradient[rep(seq_len(n), size), , drop = FALSE] +
          row_contract(result$gradient[, seq_len(q), drop = FALSE], mapped$gradient)
    }
    reduced = logsum_score(value, derivative, n)
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
  if (contracted) state$score = state$score[, seq_len(p), drop = FALSE] +
    adaptive_score(geometry, state$score[, -seq_len(p), drop = FALSE])
  if (!is.null(marginal)) {
    state$loglik = state$loglik + marginal$value
    if (differentiated) state$score = state$score + marginal$gradient
  }
  state$mode_residual = if (adaptive) geometry$mode_residual
  state$posterior = saved
  state
}

# Standalone models: fixed GH averages node partials before expanding the
# design (integrate_index_block). Adaptive GH and history use the shared kernel
# evaluator: standalone outcome rows supply values and first scores through
# its fast path, while the model's panel function supplies higher geometry
# derivatives. History bypasses the collapsed fixed-rule block to retain nodes.
integrate_standalone = function(theta, prepared, control, fixed, q, score, posterior = FALSE) {
  if (!(q > 0 && control$quadrature == 'agh') && !posterior)
    return(integrate_person_effects(theta, prepared, control, fixed, q, score,
      m = make_model_context(theta, prepared, control, fixed, indices_only = TRUE),
      block = integrate_index_block))
  integrate_person_effects(theta, prepared, control, fixed, q, score, posterior)
}

# History targets can depend on latent effects eliminated from a reduced
# likelihood. Retain the unreduced representation, including the original
# Linear person-effect maps, so nodes and derivatives describe those targets.
# GL has no supported history rule; GH/AGH weights are normalized by
# integrate_prepared after the full likelihood has been accumulated.
posterior_rule = function(theta, prepared, control, fixed = numeric()) {
  control = current_control(control)
  stop_removed_settings(control)
  if (uses_nonhistory_quadrature(control, prepared))
    stop("Posterior/history rules are unavailable for GL integration", call. = FALSE)
  control$exact_reductions = FALSE
  integrate_prepared(theta, prepared, control, fixed, score = TRUE, posterior = TRUE)$posterior
}

# Posterior modes and differentiated adaptive geometry.

row_max_abs = function(x) {
  if (ncol(x) == 1L) abs(x[, 1L]) else pmax(abs(x[, 1L]), abs(x[, 2L]))
}

panel_geometry = function(theta, prepared, control, fixed, derivatives = TRUE,
                          context = NULL, contracted = FALSE) {
  if (is.null(context)) context = make_model_context(theta, prepared, control, fixed)
  kernel = make_panel_kernel(theta, prepared, control, fixed,
    order = 2L, parameters = FALSE, latent = TRUE, context = context)
  values = make_panel_kernel(theta, prepared, control, fixed,
    order = 0L, parameters = FALSE, latent = TRUE, context = context)
  q = kernel$q; n = prepared$n_panels; p = length(theta)
  # A seed only accelerates the solve. PLN keeps its established cold-mode
  # policy: a finite inner mixture may have more than one stationary point.
  state = prepared$integration_cache$mode_state
  key = prepared$integration_cache$mode_key
  reusable = prepared$family != 'pln' && !is.null(state)
  seed = if (reusable) (state$anchor %||% state$seeds)[[key]] else NULL
  if (!is.null(seed) && (!identical(dim(seed$mode), c(n, q)) ||
      !identical(seed$parameters, names(theta)) || !identical(seed$fixed, fixed))) seed = NULL
  # Large optimizer jumps can make an old mode much worse than zero. This only
  # chooses a starting point; the solve and its acceptance criteria stay the same.
  coordinates = theta
  scales = startsWith(names(theta), 'sd_') & theta > 0
  correlations = startsWith(names(theta), 'cor_')
  coordinates[scales] = log(theta[scales])
  coordinates[correlations] = atanh(theta[correlations])
  if (!is.null(seed) && any(abs(coordinates - seed$coordinates) >
      1 + abs(seed$coordinates))) seed = NULL
  # Polish the log-concave routes to make cold and warm answers agree. Honour
  # tighter user requests; accept a representational floor only within mode_tol.
  tolerance = if (prepared$family == 'pln') control$mode_tol else min(control$mode_tol, 1e-12)
  initial = if (is.null(seed)) matrix(0, n, q) else seed$mode
  if (!is.null(state)) {
    counter = if (is.null(seed)) 'cold' else 'warm'
    state$counts[counter] = state$counts[counter] + 1L
  }
  found = tryCatch(panel_mode(kernel, values, initial, control, tolerance),
    panelselect_integration_error = function(error) {
      if (is.null(seed)) stop(error)
      state$counts['fallback'] = state$counts['fallback'] + 1L
      panel_mode(kernel, values, matrix(0, n, q), control, tolerance)
    })
  mode = found$mode; current = found$current

  full = if (derivatives) make_panel_kernel(theta, prepared, control, fixed,
    order = 3L, parameters = TRUE, latent = TRUE, context = context)$evaluate(mode, prior = TRUE) else current
  geometry = geometry_from_derivatives(full, q, p, derivatives, contracted)
  if (reusable) state$pending[[key]] = list(mode = mode, parameters = names(theta),
    fixed = fixed, coordinates = coordinates)
  c(list(mode = mode, mode_residual = found$residual, iterations = found$iterations), geometry)
}

# Safeguarded Newton solve. Values alone drive backtracking; latent derivatives
# are recomputed only at accepted points. No parameter derivatives are needed.
panel_mode = function(kernel, values, mode, control, tolerance) {
  n = nrow(mode); q = ncol(mode)
  previous_mode = mode + Inf
  current = kernel$evaluate(mode, prior = TRUE)
  done = rep(FALSE, n)
  iterations = 0L
  for (iteration in seq_len(control$mode_maxit)) {
    iterations = iteration
    if (any(!is.finite(current$value)) || any(!is.finite(current$gradient)) ||
        any(!is.finite(current$hessian)))
      integration_error("Adaptive mode has nonfinite likelihood derivatives")
    gradient = current$gradient
    scale = 1 + row_max_abs(mode)
    residual = row_max_abs(gradient)
    roundoff = row_max_abs(mode - previous_mode) <= 16*.Machine$double.eps*scale
    done = residual <= tolerance*scale | (roundoff & residual <= control$mode_tol*scale)
    if (all(done)) break
    direction = matrix(0, n, q)
    # Batched upper-Cholesky solves for the one/two person-effect dimensions.
    # Nonpositive curvature keeps the gradient direction, as in the scalar solve.
    active = which(!done)
    step = gradient[active, , drop = FALSE]
    H11 = -current$hessian[active, 1L, 1L]
    positive = H11 > 0
    r11 = sqrt(H11[positive])
    if (q == 1L) step[positive, 1L] = (step[positive, 1L]/r11)/r11
    else {
      eligible = which(positive)
      r12 = -current$hessian[active[eligible], 1L, 2L]/r11
      remaining = -current$hessian[active[eligible], 2L, 2L]-r12*r12
      keep = which(remaining > 0)
      eligible = eligible[keep]; r11 = r11[keep]; r12 = r12[keep]
      r22 = sqrt(remaining[keep])
      first = step[eligible, 1L]/r11
      second = ((step[eligible, 2L]-r12*first)/r22)/r22
      step[eligible, 1L] = (first-r12*second)/r11
      step[eligible, 2L] = second
    }
    direction[active, ] = step/pmax(1, sqrt(rowSums(step*step))/5)
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
    previous_mode = mode
    mode = proposal
    current = kernel$evaluate(mode, prior = TRUE)
  }
  residual = row_max_abs(current$gradient)
  if (any(!is.finite(residual)) ||
      any(residual > 10*control$mode_tol*(1 + row_max_abs(mode))))
    integration_error("Adaptive mode did not satisfy its gradient tolerance")

  list(mode = mode, current = current, residual = residual, iterations = iterations)
}

# Upper-Cholesky map C = chol(H)^(-1), H = -h_zz, batched over panels.
# Forward derivatives are needed for posterior nodes. Fitting instead retains
# the mixed derivatives at the mode and contracts them after reducing nodes.
geometry_from_derivatives = function(full, q, p, derivatives, contracted = FALSE) {
  n = length(full$value)
  H11 = -full$hessian[, 1L, 1L]
  if (any(!is.finite(H11) | H11 <= 0))
    integration_error('Adaptive mode curvature is not positive definite')
  r11 = sqrt(H11)
  transform = array(0, c(n, q, q)); transform[, 1L, 1L] = 1/r11
  log_jacobian = -log(r11)
  if (q == 2L) {
    r12 = -full$hessian[, 1L, 2L]/r11
    H22_conditional = -full$hessian[, 2L, 2L] - r12*r12
    if (any(!is.finite(H22_conditional) | H22_conditional <= 0))
      integration_error('Adaptive mode curvature is not positive definite')
    r22 = sqrt(H22_conditional)
    c12 = -r12/(r11*r22)
    transform[, 1L, 2L] = c12; transform[, 2L, 2L] = 1/r22
    log_jacobian = log_jacobian - log(r22)
  }
  out = list(transform = transform, log_jacobian = log_jacobian)
  if (!derivatives) return(out)
  if (contracted) return(c(out, list(hessian = full$hessian, third = full$third)))
  dm = array(0, c(n, q, p)); dC = array(0, c(n, q, q, p))
  columns = q + seq_len(p)
  if (q == 1L) {
    dm1 = matrix(full$hessian[, 1L, columns], n, p)/H11
    dH = -matrix(full$third[, 1L, 1L, columns], n, p) - dm1*full$third[, 1L, 1L, 1L]
    dm[, 1L, ] = dm1
    dC[, 1L, 1L, ] = -.5*dH/(H11*r11)
    dlog_jacobian = -.5*dH/H11
  } else {
    g1 = matrix(full$hessian[, 1L, columns], n, p)
    g2 = matrix(full$hessian[, 2L, columns], n, p)
    dm2 = (g2 - r12*g1/r11)/(r22*r22)
    dm1 = g1/H11 - r12*dm2/r11
    dm[, 1L, ] = dm1; dm[, 2L, ] = dm2
    dH = function(a, b) -matrix(full$third[, a, b, columns], n, p) -
      dm1*full$third[, a, b, 1L] - dm2*full$third[, a, b, 2L]
    dr11 = dH(1L, 1L)/(2*r11)
    dr12 = (dH(1L, 2L) - r12*dr11)/r11
    dr22 = (dH(2L, 2L) - 2*r12*dr12)/(2*r22)
    dC[, 1L, 1L, ] = -dr11/H11
    dC[, 2L, 2L, ] = -dr22/(r22*r22)
    dC[, 1L, 2L, ] = -dr12/(r11*r22) - c12*(dr11/r11 + dr22/r22)
    dlog_jacobian = -dr11/r11 - dr22/r22
  }
  c(out, list(mode_gradient = dm, transform_gradient = dC,
    jacobian_gradient = dlog_jacobian))
}

# Complete moving-node/determinant correction after quadrature reduction.
# moments contain E[h_z] followed by E[h_z1*x1], and (q=2) E[h_z1*x2], E[h_z2*x2].
# Only upper-triangular C entries move. Curvature weights below use the unique
# H11, H12, H22 coordinates: H12 represents both symmetric entries, once.
adaptive_score = function(geometry, moments) {
  n = nrow(moments); q = ncol(geometry$mode)
  hessian = geometry$hessian; third = geometry$third
  p = dim(hessian)[3L] - q; columns = q + seq_len(p)
  C = geometry$transform; c11 = C[, 1L, 1L]
  if (q == 1L) {
    S11 = -.5*c11*c11*(1 + moments[, 2L]*c11)
    v1 = (moments[, 1L] - S11*third[, 1L, 1L, 1L])*c11*c11
    return(v1*matrix(hessian[, 1L, columns], n, p) -
      S11*matrix(third[, 1L, 1L, columns], n, p))
  }
  c12 = C[, 1L, 2L]; c22 = C[, 2L, 2L]
  r11 = 1/c11; r22 = 1/c22; r12 = -c12*r11*r22
  # Reverse the same triangular formulas used in geometry_from_derivatives().
  # Adding C^(-T) supplies the log-determinant derivative.
  G11 = moments[, 3L] + r11; G12 = moments[, 4L]; G22 = moments[, 5L] + r22
  dr11 = -(G11*c11 + G12*c12)/r11
  dr12 = -G12/(r11*r22)
  dr22 = -(G22*c22 + G12*c12)/r22
  S22 = dr22/(2*r22)
  S12 = (dr12 - 2*r12*S22)/r11
  S11 = (dr11 - r12*S12)/(2*r11)
  b1 = moments[, 1L] - S11*third[, 1L, 1L, 1L] -
    S12*third[, 1L, 2L, 1L] - S22*third[, 2L, 2L, 1L]
  b2 = moments[, 2L] - S11*third[, 1L, 1L, 2L] -
    S12*third[, 1L, 2L, 2L] - S22*third[, 2L, 2L, 2L]
  v2 = (b2 - r12*b1/r11)/(r22*r22)
  v1 = b1/(r11*r11) - r12*v2/r11
  v1*matrix(hessian[, 1L, columns], n, p) + v2*matrix(hessian[, 2L, columns], n, p) -
    S11*matrix(third[, 1L, 1L, columns], n, p) -
    S12*matrix(third[, 1L, 2L, columns], n, p) - S22*matrix(third[, 2L, 2L, columns], n, p)
}
