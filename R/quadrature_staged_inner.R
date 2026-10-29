# Staged GL for the PLN observation effect. Each active PLN
# outcome row has its own positive observation mesh, fixed within an optimization
# stage (only selected rows for selection models). Model files supply the
# observation kernel and local partials; this file evaluates it
# on canonical and reference rules, records inner errors and tails, and refines
# a row's mesh from its probes between stages (prepare_staged_integration).
# Inner GL rules include the normal density in their weights. Outer GL mesh
# rules instead request it from the panel evaluator; do not apply either measure twice.

# 1. Observation rules and diagnostic records ----------------------------

staged_inner_rule = function(mesh, reference = FALSE) {
  rule = staged_cell_rule(mesh, reference = reference)
  # These are weights for the normal expectation, unlike the outer rule where
  # the model kernel explicitly includes the normal density.
  rule$log_weights = rule$log_weights+dnorm(rule$nodes[, 1L], log = TRUE)
  rule
}

staged_inner_tail = function(A, B, y, gamma, tau, mesh, value, gradient, parameters) {
  bound = mesh$bound; logs = staged_tail_logs(bound, 1L)
  bounds = staged_count_bounds(y); selected = !is.null(A); c_tau = sqrt(1-tau^2)
  mass = exp(log(bounds$mass)+logs[['mass']]-value)
  first = exp(logs[['first']]-value)
  numerator = c(if (selected) bounds$mass*dnorm(0)/c_tau*exp(logs[['mass']]-value),
    bounds$first*exp(logs[['mass']]-value), if (parameters) bounds$first*first,
    if (selected && parameters) bounds$mass*(dnorm(0)/c_tau*first+
      abs(tau)*dnorm(1)/c_tau^2*exp(logs[['mass']]-value)))
  list(relative = mass, score = numerator+abs(drop(gradient))*mass)
}

staged_log_difference = function(a, b) {
  high = pmax(a, b)
  result = high+log(-expm1(-abs(a-b)))
  result[a == b] = -Inf
  result
}

staged_inner_probes = function(recorder, ids = seq_along(recorder$inner_probes)) {
  probes = list()
  for (id in ids) {
    value = recorder$inner_probes[[id]]; derivative = recorder$inner_score_probes[[id]]
    if (!is.null(value)) probes[[length(probes)+1L]] = value
    if (!is.null(derivative) && !identical(value, derivative))
      probes[[length(probes)+1L]] = derivative
  }
  probes
}

# Log error estimates for a batch of PLN observation rows: the empirical
# canonical/reference discrepancies plus analytic omitted-tail bounds. The
# finite-rule discrepancies do not bound the continuous quadrature remainder.
# Dimensions: rows x 1 probability logs; rows x derivatives first logs.
# Batching shares only these diagnostic reductions; nodes, model kernels,
# finite-objective derivatives and accuracy budgets are those of each row.
staged_inner_error_logs_batch = function(canonical, reference, y, tau, mesh,
                                          selected, parameters, order) {
  n = length(canonical$value)
  logs = staged_tail_logs(mesh$bound, 1L)
  distinct = unique(y); index = match(y, distinct)
  bounds = staged_count_bounds(distinct)
  mass = bounds$mass[index]; first_bound = bounds$first[index]
  probability = as.vector(logsum_weights(c(
    staged_log_difference(canonical$value, reference$value),
    log(mass)+logs[['mass']]), n)$value)
  if (!order) return(list(probability = probability,
    first = matrix(numeric(), n, 0L)))

  high = pmax(canonical$value, reference$value)
  # Matrix recycling is row-fast, preserving each scalar query's arithmetic.
  first = high+log(abs(exp(canonical$value-high)*canonical$gradient-
    exp(reference$value-high)*reference$gradient))
  d = ncol(first); tails = matrix(0, n, d); j = 1L
  c_tau = sqrt(1-tau^2)
  if (selected) {
    tails[, j] = log(mass*dnorm(0)/c_tau)+logs[['mass']]
    j = j+1L
  }
  tails[, j] = log(first_bound)+logs[['mass']]; j = j+1L
  if (parameters) {
    tails[, j] = log(first_bound)+logs[['first']]; j = j+1L
  }
  if (selected && parameters) tails[, j] = log(mass)+
    as.vector(logsum_weights(c(log(dnorm(0)/c_tau)+logs[['first']],
      log(abs(tau)*dnorm(1)/c_tau^2)+logs[['mass']]), 1L)$value)
  first = matrix(logsum_weights(c(as.vector(first), as.vector(tails)), n*d)$value,
    n, d)
  list(probability = probability, first = first)
}

staged_record_inner_batch = function(recorder, rule, id, A, B, y, gamma, tau,
                                     log_factorial, canonical, reference, order,
                                     parameters, batch_nodes = 4096L) {
  n = length(B); if (!n) return(invisible(NULL))
  selected = !is.null(A)
  d = if (order) ncol(canonical$gradient) else 0L
  # No N-by-d-by-p tensor. At most two N-by-d reduction buffers per natural
  # parameter; cap independently of the caller's kernel-query batch size.
  size = max(1L, min(as.integer(batch_nodes), floor(buffer_bytes/(8*32*max(1L, d)))))
  names_local = c(if (selected) 'a', 'b', if (parameters) 'error',
    if (selected && parameters) 'tau')
  base = rule$projection; local = match(id, rule$row_ids)
  p = length(rule$theta); names_theta = names(rule$theta)
  target = match(names_theta, recorder$parameter_names)
  free = !names_theta %in% recorder$fixed_names
  b = match('b', names_local); a = match('a', names_local)
  lambda = match(base$lambda_name, names_theta)
  delta = match('sd_re_selection', names_theta)
  rho = match('cor_re', names_theta)
  fixed_map = matrix(0, d, p)
  if (order) {
    if (length(base$beta)) fixed_map[b, base$beta] = abs(base$x[local, ])
    if (selected && length(base$alpha)) fixed_map[a, base$alpha] = abs(base$w[local, ])
    if (parameters) {
      k = match('sd_error', names_theta)
      if (!is.na(k)) fixed_map[match('error', names_local), k] = 1
      k = match('cor_error', names_theta)
      if (!is.na(k) && selected) fixed_map[match('tau', names_local), k] = 1
    }
  }
  probe = function(k) list(row = id, A = if (selected) A[k] else NULL,
    B = B[k], y = y[k], log_factorial = log_factorial[k])

  for (start in seq.int(1L, n, by = size)) {
    take = seq.int(start, min(n, start+size-1L)); count = length(take)
    canonical_batch = list(value = canonical$value[take],
      gradient = if (order) canonical$gradient[take, , drop = FALSE] else NULL)
    reference_batch = list(value = reference$value[take],
      gradient = if (order) reference$gradient[take, , drop = FALSE] else NULL)
    error = staged_inner_error_logs_batch(canonical_batch, reference_batch, y[take], tau,
      rule$meshes[[id]], selected, parameters, order)
    if (anyNA(error$probability)) {
      bad = which(is.na(error$probability))
      probes = lapply(take[bad], probe)
      staged_mesh_required('PLN inner likelihood discrepancy is not finite', rule$theta,
        list(row = id, order = order, parameters = parameters,
          # current_value keeps its recorded name in failure diagnostics.
          current_value = canonical_batch$value[bad], reference_value = reference_batch$value[bad],
          y = y[take[bad]], A = if (selected) A[take[bad]] else NULL,
          B = B[take[bad]], log_factorial = log_factorial[take[bad]],
          mesh = rule$meshes[[id]][c('method','bound')]),
        list(list(theta = rule$theta, inner_probes = probes)))
    }
    maximum = max(error$probability)
    # Scalar >= selects the last equal maximum, including an all -Inf batch.
    if (maximum >= recorder$inner_error[id]) {
      k = tail(which(error$probability == maximum), 1L)
      recorder$inner_probes[[id]] = probe(take[k])
    }
    recorder$inner_error[id] = max(recorder$inner_error[id], maximum)
    if (!order || !ncol(error$first)) next

    V = if (base$lambda > 0) (B[take]-base$b[local])/base$lambda else rep(0, count)
    U = if (selected && base$delta > 0)
      (A[take]-base$a[local])/base$delta else rep(0, count)
    candidate = rep(-Inf, count)
    for (k in seq_len(p)) {
      map = matrix(rep(fixed_map[, k], each = count), count, d)
      if (!is.na(lambda) && k == lambda) map[, b] = abs(V)
      if (selected && !is.na(delta) && k == delta) map[, a] = abs(U)
      if (!is.na(rho) && k == rho && base$delta > 0 && base$lambda > 0) {
        if (selected) {
          Z = (V-base$rho*U)/sqrt(1-base$rho^2)
          map[, b] = abs(base$lambda*(U-base$rho/sqrt(1-base$rho^2)*Z))
        } else map[, b] = base$lambda*recorder$outer_radius*
          (1+abs(base$rho)/sqrt(1-base$rho^2))
      }
      error_k = as.vector(logsum_weights(as.vector(error$first+log(map)), count)$value)
      recorder$inner_score_error[id, target[k]] = max(
        recorder$inner_score_error[id, target[k]], error_k)
      if (free[k]) candidate = pmax(candidate, error_k)
    }
    maximum = max(candidate)
    if (maximum > -Inf && maximum >= recorder$inner_score_probe_error[id]) {
      k = tail(which(candidate == maximum), 1L)
      recorder$inner_score_probes[[id]] = probe(take[k])
      recorder$inner_score_probe_error[id] = maximum
    }
  }
  invisible(NULL)
}

staged_inner_tail_relative_batch = function(y, mesh, value) {
  distinct = unique(y); bounds = staged_count_bounds(distinct)
  exp(log(bounds$mass[match(y, distinct)])+
    staged_tail_logs(mesh$bound, 1L)[['mass']]-value)
}

integrate_pln_inner_staged = function(A, B, y, gamma, tau, order, log_factorial,
                                     rule, kernel, parameters = TRUE) {
  selected = !is.null(A); n = length(B)
  coordinates = c(if (selected) 'a', 'b', if (parameters) 'error',
    if (selected && parameters) 'tau')
  r = if (selected) 2L else 1L
  out = local_block(numeric(n), coordinates, r, order)
  if (!n) return(out)
  paired = isTRUE(rule[['paired']])
  if (paired) {
    if (is.null(rule$recorder) || isTRUE(rule[['reference']]))
      integration_error('Paired PLN inner evaluation requires canonical recording')
    reference_out = local_block(numeric(n), coordinates, r, order)
  }
  if (is.null(rule$meshes) || !length(rule$row_ids))
    staged_mesh_required('PLN observation mesh is unavailable', rule$theta)
  ids = rep_len(rule$row_ids, n)
  recorder = rule$recorder
  budget = rule$budget
  if (!is.null(budget)) on.exit(budget$end(check_limit = FALSE), add = TRUE)
  inner_clock = function() {
    rule$clock()
    if (!is.null(budget)) budget$check()
    invisible(NULL)
  }
  for (id in unique(ids)) {
    if (!is.null(budget)) budget$begin(id)
    take = which(ids == id); mesh = rule$meshes[[id]]
    if (is.null(mesh)) staged_mesh_required('PLN row has no retained mesh', rule$theta,
      list(row = id, A = if (selected) A[take[1L]] else NULL, B = B[take[1L]]))
    staged_resource_check(mesh, rule$options, inner_clock)
    # All queries for this row (including new outer GH mode nodes) share these
    # fixed abscissae. No query can lazily change a row's mesh.
    weights = staged_inner_rule(mesh, isTRUE(rule[['reference']]))
    reference_rule = if (!is.null(recorder) && !isTRUE(rule[['reference']]))
      staged_inner_rule(mesh, TRUE) else NULL
    if (!is.null(reference_rule)) {
      recorder$inner_mass[id] = max(1, recorder$inner_mass[id],
        exp(logsum_weights(weights$log_weights, 1L)$value),
        exp(logsum_weights(reference_rule$log_weights, 1L)$value))
      recorder$inner_first[id] = max(recorder$inner_first[id],
        sum(exp(weights$log_weights)*abs(weights$nodes[, 1L])),
        sum(exp(reference_rule$log_weights)*abs(reference_rule$nodes[, 1L])))
    }
    nodes = length(weights$log_weights); d = length(coordinates)
    entries = 1+d*(order > 0L)+r*d*(order > 1L)+r*r*d*(order > 2L)
    size = max(1L, min(rule$options$batch_nodes,
      floor(buffer_bytes/(8*6*nodes*entries))))
    for (first in seq.int(1L, length(take), by = size)) {
      rows = take[seq.int(first, min(length(take), first+size-1L))]; count = length(rows)
      canonical_part = observation_mixture(kernel, A, B, y, gamma, tau, rows,
        weights$nodes[, 1L], weights$log_weights, order, log_factorial, parameters, r)
      if (!is.null(recorder)) recorder$inner_work[id] = recorder$inner_work[id]+count*nodes
      out$value[rows] = canonical_part$value
      if (order) out$gradient[rows, ] = canonical_part$gradient
      if (order >= 2L) out$hessian[rows, , ] = canonical_part$hessian
      if (order >= 3L) out$third[rows, , , ] = canonical_part$third
      if (!is.null(recorder) && !isTRUE(rule[['reference']])) {
        tail = staged_inner_tail_relative_batch(y[rows], mesh, canonical_part$value)
        recorder$inner_tail[id] = max(recorder$inner_tail[id], tail)
        reference_part = observation_mixture(kernel, A, B, y, gamma, tau, rows,
          reference_rule$nodes[, 1L], reference_rule$log_weights, order, log_factorial,
          parameters, r)
        recorder$inner_work[id] = recorder$inner_work[id]+count*length(reference_rule$log_weights)
        staged_record_inner_batch(recorder, rule, id,
          if (selected) A[rows] else NULL, B[rows], y[rows], gamma, tau,
          log_factorial[rows], canonical_part, reference_part, order, parameters)
        if (paired) {
          reference_out$value[rows] = reference_part$value
          if (order) reference_out$gradient[rows, ] = reference_part$gradient
          if (order >= 2L) reference_out$hessian[rows, , ] = reference_part$hessian
          if (order >= 3L) reference_out$third[rows, , , ] = reference_part$third
        }
      }
      inner_clock()
    }
    if (!is.null(budget)) budget$end()
  }
  if (any(!is.finite(out$value)) || (order && any(!is.finite(out$gradient))))
    staged_mesh_required('PLN observation likelihood or derivative is not finite', rule$theta,
      list(rows = unique(ids), inner_probes = if (!is.null(recorder))
        staged_inner_probes(recorder, unique(ids)) else NULL))
  if (paired) {
    if (any(!is.finite(reference_out$value)) ||
        (order && any(!is.finite(reference_out$gradient))))
      staged_mesh_required('Paired PLN reference likelihood or derivative is not finite', rule$theta,
        list(rows = unique(ids), inner_probes = staged_inner_probes(recorder, unique(ids))))
    out$reference_block = reference_out
  }
  out
}

# 2. Refine observation meshes between stages ----------------------------

staged_inner_probe_value = function(theta, prepared, control, fixed, mesh, probe, clock) {
  selected = !is.null(probe$A)
  kernel = if (selected) probit_re_pln_re_inner_kernel else pln_re_inner_kernel
  d = if (selected) 4L else 2L
  canonical_rule = staged_inner_rule(mesh)
  reference_rule = staged_inner_rule(mesh, TRUE)
  evaluate = function(rule) {
    n = nrow(rule$nodes)
    term = kernel(if (selected) rep(probe$A, n) else NULL, rep(probe$B, n),
      rep(probe$y, n), theta[['sd_error']], if (selected) theta[['cor_error']] else 0,
      rule$nodes[, 1L], rep(0, n), 1L, rep(probe$log_factorial, n), TRUE)
    staged_reduce_cells(term, rule$log_weights, rule$cells, d)
  }
  clock(); canonical = evaluate(canonical_rule); reference = evaluate(reference_rule)
  canonical_sum = logsum_score(canonical[, 1L], canonical[, -1L, drop = FALSE], 1L)[1L, ]
  reference_sum = logsum_score(reference[, 1L], reference[, -1L, drop = FALSE], 1L)[1L, ]
  if (!is.finite(canonical_sum[1L]) || !is.finite(reference_sum[1L]))
    return(list(finite = FALSE, cell_error = rep(Inf, nrow(mesh$lower))))
  w = exp(canonical[, 1L]-canonical_sum[1L]); wr = exp(reference[, 1L]-canonical_sum[1L])
  cell_error = abs(w-wr)
  first_error = staged_score_discrepancy(canonical[, -1L, drop = FALSE],
    reference[, -1L, drop = FALSE], w, wr, canonical_sum[-1L], reference_sum[-1L])
  tail = staged_inner_tail(probe$A, probe$B, probe$y, theta[['sd_error']],
    if (selected) theta[['cor_error']] else 0, mesh, canonical_sum[1L], canonical_sum[-1L], TRUE)
  list(finite = TRUE, loglik = canonical_sum[1L], score = canonical_sum[-1L],
    loglik_error = max(abs(canonical_sum[1L]-reference_sum[1L]), sum(cell_error)),
    score_error = first_error+tail$score, tail = tail$relative,
    tail_score = tail$score,
    cell_error = cell_error+rowSums(abs(canonical[, -1L, drop = FALSE]*w-
      reference[, -1L, drop = FALSE]*wr)))
}

staged_prepare_inner = function(theta, prepared, control, fixed, mesh, probe, clock,
                                tightening = 1, absolute_target = NULL, tail_absolute_target = NULL) {
  options = staged_options(control, TRUE)
  budget = control$.staged_inner_budget
  if (!is.null(budget)) {
    on.exit(budget$end(check_limit = FALSE), add = TRUE)
    budget$begin(probe$row)
  }
  inner_clock = function() {
    clock()
    if (!is.null(budget)) budget$check()
    invisible(NULL)
  }
  for (iteration in seq_len(staged_refinement_passes)) {
    staged_resource_check(mesh, options, inner_clock, staged_options_name(control, TRUE))
    value = staged_inner_probe_value(theta, prepared, control, fixed, mesh, probe, inner_clock)
    inner_clock()
    if (!value$finite)
      staged_mesh_required('PLN observation mesh has no finite positive mass', theta,
        list(inner_probes = list(probe)))
    tail_target = options$tail_tol
    if (!is.null(tail_absolute_target)) tail_target = min(tail_target,
      exp(tail_absolute_target-value$loglik))
    good = value$loglik_error <= options$panel_loglik_tol/tightening &&
      all(value$score_error <= options$panel_score_tol/tightening*(1+abs(value$score))) &&
      value$tail <= tail_target
    if (!is.null(absolute_target)) good = good &&
      log(value$loglik_error)+value$loglik <= absolute_target &&
      all(log(value$score_error)+value$loglik <= absolute_target)
    if (good) {
      if (!is.null(budget)) budget$end()
      return(mesh)
    }
    if (!is.null(absolute_target)) tail_target = min(tail_target,
      exp(absolute_target-value$loglik)/(1+max(abs(value$score))))
    # Subdivision cannot shrink an omitted derivative tail. Its normalized
    # and absolute-numerator budgets are the same ones used above to accept
    # this row. The preparer checks all local a/b/error/tau
    # derivatives, including fixed occasion parameters, for latent/mixed
    # calculations. Global fixed-coordinate masking belongs to panel gates.
    tail_score_budget = options$panel_score_tol/tightening*(1+abs(value$score))
    if (!is.null(absolute_target)) tail_score_budget = pmin(tail_score_budget,
      exp(absolute_target-value$loglik))
    score_tail_bad = any(!is.finite(value$tail_score) |
      value$tail_score > tail_score_budget)
    if (value$tail > tail_target || score_tail_bad)
      staged_bound_failure(if (score_tail_bad)
          'PLN observation omitted-score-tail envelope unresolved at maximum bound' else
          'PLN observation omitted-tail envelope unresolved at maximum bound',
        options$bound, staged_options_name(control, TRUE, bound = TRUE), theta,
        list(inner_probes = list(probe),
          component = 'inner', bound = mesh$bound, tail = value$tail,
          score_tail = value$tail_score, score_budget = tail_score_budget,
          score_coordinates_checked = 'all original local a/b/error/tau coordinates'))
    mesh = staged_refine_cells(mesh, staged_select_cells(value$cell_error),
      rep(1L, length(staged_select_cells(value$cell_error))))
  }
  staged_mesh_required(paste0('PLN observation mesh refinement limit exhausted (',
    staged_refinement_passes, ' refinements of row ', probe$row, ')'), theta,
    list(inner_probes = list(probe)))
}

# Tighten each recorded inner probe to an equal share of an absolute panel
# target; companion factors and the outer measure convert it to one row's error.
staged_refine_inner_probes = function(inner, probes, source, loglik, target, tail_target,
                                      theta, prepared, control, fixed, clock, tightening) {
  points = Filter(Negate(is.null), probes)
  for (probe in points) {
    companion = source$inner_companion_log[as.character(probe$row)] %||% 0
    measure = max(1, source$inner_outer_mass %||% 1)
    inner[[probe$row]] = staged_prepare_inner(theta, prepared, control, fixed,
      inner[[probe$row]], probe, clock, tightening,
      absolute_target = loglik+log(target/(8*max(1, length(points))))-companion-log(measure),
      tail_absolute_target = loglik+log(tail_target/(8*max(1, length(points))))-
        companion-log(measure))
  }
  inner
}

# Seed or reuse each active row's observation mesh, then prepare it at theta.
staged_seed_inner = function(mesh, theta, prepared, control, fixed, context, previous,
                             clock, tightening) {
  active = if (prepared$has_selection) context$data$selected$rows else context$data$active$rows
  for (row in active) {
    current = if (!is.null(previous$inner[[row]])) previous$inner[[row]] else
      staged_seed_mesh(context, prepared, theta, control, fixed, inner = TRUE, row = row)
    current$method = control$observation_quadrature
    inner_options = staged_options(control, TRUE)
    current$rule = staged_rule(control$observation_nodes, inner_options$reference_nodes, 1L)
    # A changed user bound is honored on a diagnostic mesh; a tighter bound
    # never discards a previously retained interior merely for speed.
    if (current$bound < inner_options$bound) current = staged_extend_cells(current, inner_options$bound)
    probe = list(row = row, A = if (prepared$has_selection && !context$tau_zero) context$a[row] else NULL,
      B = context$b[row], y = prepared$y[row], log_factorial = context$data$log_factorial[row])
    mesh$inner[[row]] = staged_prepare_inner(theta, prepared, control, fixed, current, probe,
      clock, tightening)
  }
  mesh
}
