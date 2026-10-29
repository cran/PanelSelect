# Selected panel PLN with occasion correlation tau fixed at zero. This route
# requires GL at both person and observation levels, positive effect scales,
# and the remaining eligibility checks below. A free tau evaluated at zero
# does not qualify, because its coupled derivative must remain available.
#
# Frozen positive rules in physical (u,v) coordinates reuse each outcome
# integral across selection nodes: the selection factor depends on u, the
# outcome factor on v, and their normal density carries rho. R supplies the
# model/prior scores; generic compiled reducers only sum positive measures.

staged_physical_pln_eligible = function(theta, prepared, control, fixed) {
  prepared$family == 'pln' && prepared$has_selection &&
    !identical(prepared$structure, 'cross_sectional') &&
    !identical(control$.physical_pln, FALSE) &&
    (control$observation_quadrature %||% 'gh') == 'gl' && uses_gl_quadrature(control) &&
    is_fixed_zero('cor_error', theta, fixed) && isTRUE(theta[['cor_error']] == 0) &&
    all(theta[c('sd_re_selection','sd_re_outcome','sd_error')] > 0) &&
    all(prepared$active %||% rep(TRUE, prepared$n_rows))
}

staged_physical_node_count = function(mesh) {
  u = nrow(mesh$u$lower); v = nrow(mesh$v$lower)
  u*v*(length(mesh$u$rule$nodes)*length(mesh$v$rule$nodes)+
    length(mesh$u$rule[['reference']]$nodes)*length(mesh$v$rule[['reference']]$nodes))
}

staged_physical_refine = function(mesh, value) {
  pairs = as.matrix(expand.grid(seq_len(nrow(mesh$u$lower)), seq_len(nrow(mesh$v$lower))))
  cells = staged_select_cells(value$cell_refinement_error)
  dimension = max.col(value$dimension_error[cells, , drop = FALSE], ties.method = 'first')
  for (j in 1:2) {
    take = unique(pairs[cells[dimension == j], j])
    if (length(take)) {
      axis = if (j == 1L) 'u' else 'v'
      mesh[[axis]] = staged_refine_cells(mesh[[axis]], take, rep(1L, length(take)))
    }
  }
  mesh
}

staged_physical_prior = function(u, v, rho, factors = NULL,
                               score = TRUE, absolute_cross = FALSE) {
  q = 1-rho^2
  a = u$log_weights-u$nodes[, 1L]^2/(2*q)
  b = v$log_weights-v$nodes[, 1L]^2/(2*q)
  if (!is.null(factors)) {a = a+factors$selection; b = b+factors$outcome}
  U = u$nodes[, 1L]; V = v$nodes[, 1L]
  out = bilinear_reduce_cpp(a, b, U, V, rho/q, score)
  out$value = out$value-log(2*pi)-.5*log(q)
  if (score) {
    out$first_u = sum(out$mu*abs(U)); out$first_v = sum(out$mv*abs(V))
    out$second_u = sum(out$mu*U^2); out$second_v = sum(out$mv*V^2)
    out$cross = sum(U*out$weighted_v)
    out$rho_score = rho/q+((1+rho^2)*out$cross-
      rho*(out$second_u+out$second_v))/q^2
    # Actual finite normal measure; absolute density-score envelope uses its
    # measured second moments, never exact-normal moments or assumed mass1.
    out$absolute_cross_bound = .5*(out$second_u+out$second_v)
    if (absolute_cross) {
      out$absolute_cross = 0
      # Four disjoint nonzero quadrants reuse the generic streaming reducer.
      # Zero coordinates contribute zero to |UV| and are omitted only here.
      for (su in c(-1, 1)) for (sv in c(-1, 1)) {
        iu = which(sign(U) == su); iv = which(sign(V) == sv)
        if (!length(iu) || !length(iv)) next
        uq = list(nodes = u$nodes[iu, , drop = FALSE], log_weights = u$log_weights[iu])
        vq = list(nodes = v$nodes[iv, , drop = FALSE], log_weights = v$log_weights[iv])
        fq = if (is.null(factors)) NULL else list(selection = factors$selection[iu],
          outcome = factors$outcome[iv])
        quadrant = staged_physical_prior(uq, vq, rho, fq, TRUE, FALSE)
        out$absolute_cross = out$absolute_cross+
          exp(quadrant$value-out$value)*abs(quadrant$cross)
      }
      out$absolute_cross_bound = out$absolute_cross
    }
    out$rho_first_bound = abs(rho)/q+(1+rho^2)*out$absolute_cross_bound/q^2+
      abs(rho)*(out$second_u+out$second_v)/q^2
  }
  out
}

staged_physical_pln_diagnostics = function(theta, prepared, fixed, mesh, recorder,
                                     prior, loglik, score, derivatives = TRUE) {
  p = length(theta); selected = which(prepared$d == 1); y = prepared$y[selected]
  ids = (prepared$.staged_rows %||% seq_len(prepared$n_rows))[selected]
  bounds = staged_count_bounds(y); log_maximum = log(bounds$mass)
  log_C = log_maximum+log(recorder$inner_mass[ids])
  companion = sum(log_C)-log_C
  measure = exp(prior$value)
  probability = measure*sum(exp(recorder$inner_error[ids]+companion-loglik))
  if (!derivatives) {
    log_factor = sum(log_maximum)-loglik
    log_mass = as.numeric(logsum_weights(c(staged_tail_logs(mesh$u$bound, 1L)[['mass']],
      staged_tail_logs(mesh$v$bound, 1L)[['mass']]), 1L)$value)
    inner_mass = vapply(ids, function(id) staged_tail_logs(mesh$inner[[id]]$bound, 1L)[['mass']], numeric(1L))
    outer_log_tail = log_factor+log_mass
    inner_log_tail = log_factor+as.numeric(logsum_weights(inner_mass, 1L)$value)
    outer_tail = exp(outer_log_tail); inner_tail = exp(inner_log_tail)
    zero = setNames(numeric(p), names(theta))
    return(list(nested_error = probability, nested_score_error = zero,
      inner_companion_log = setNames(companion, as.character(ids)), inner_outer_mass = measure,
      inner_probes = staged_inner_probes(recorder, ids),
      outer_tail = outer_tail, outer_log_tail = outer_log_tail, outer_tail_score = zero,
      inner_tail = inner_tail, inner_log_tail = inner_log_tail, inner_tail_score = zero,
      tail = outer_tail+inner_tail, physical_prior = list(mass = measure,
        score_moments_evaluated = FALSE),
      evidence = 'empirical finite-rule probability discrepancies plus analytic omitted true-model tails'))
  }
  first = measure*colSums(exp(recorder$inner_score_error[ids, , drop = FALSE]+
    companion-loglik))
  kernel = setNames(numeric(p), names(theta))
  s_names = paste0('selection.', colnames(prepared$w))
  y_names = paste0('outcome.', colnames(prepared$x))
  kernel[s_names] = dnorm(0)*colSums(abs(prepared$w))
  kernel['sd_re_selection'] = dnorm(0)*nrow(prepared$w)*prior$first_u
  local = bounds$first*recorder$inner_mass[ids]/exp(log_C)
  kernel[y_names] = colSums(abs(prepared$x[selected, , drop = FALSE])*local)
  kernel['sd_re_outcome'] = sum(local)*prior$first_v
  kernel['sd_error'] = sum(bounds$first*recorder$inner_first[ids]/exp(log_C))
  kernel['cor_re'] = prior$rho_first_bound
  nested_score = first+probability*(kernel+abs(score))

  # Strict TRUE-MODEL omitted tails; finite component masses do not stand in
  # for true inside probabilities. All output ratios use the calculated L.
  Bu = mesh$u$bound; Bv = mesh$v$bound; rho = theta[['cor_re']]; q = 1-rho^2
  ls = function(x) as.numeric(logsum_weights(x, 1L)$value)
  ulogs = staged_tail_logs(Bu, 1L); vlogs = staged_tail_logs(Bv, 1L)
  log_mass = ls(c(ulogs[['mass']], vlogs[['mass']]))
  log_U = ls(c(ulogs[['first']], log(abs(rho))+vlogs[['first']],
    .5*log(q)+.5*log(2/pi)+vlogs[['mass']]))
  log_V = ls(c(vlogs[['first']], log(abs(rho))+ulogs[['first']],
    .5*log(q)+.5*log(2/pi)+ulogs[['mass']]))
  log_second_u = ls(c(log(Bu)+ulogs[['first']], ulogs[['mass']]))
  log_second_v = ls(c(log(Bv)+vlogs[['first']], vlogs[['mass']]))
  log_second = ls(c(log(1+rho^2)+ls(c(log_second_u, log_second_v)), log(q)+log_mass))
  log_rho = ls(c(log(abs(rho))-log(q)+log_mass,
    2*log(1+abs(rho))-log(2)-2*log(q)+log_second))
  log_factor = sum(log_maximum)-loglik
  outer_log_tail = log_factor+log_mass
  outer_tail = exp(outer_log_tail)
  outer_first = setNames(numeric(p), names(theta))
  outer_first[s_names] = exp(log_factor+log(dnorm(0)*colSums(abs(prepared$w)))+log_mass)
  outer_first['sd_re_selection'] = exp(log_factor+log(dnorm(0)*nrow(prepared$w))+log_U)
  local = bounds$first/bounds$mass
  y_coefficient = colSums(abs(prepared$x[selected, , drop = FALSE])*local)
  outer_first[y_names] = exp(log_factor+log(y_coefficient)+log_mass)
  outer_first['sd_re_outcome'] = exp(log_factor+log(sum(local))+log_V)
  outer_first['sd_error'] = exp(log_factor+log(sum(local))+.5*log(2/pi)+log_mass)
  outer_first['cor_re'] = exp(log_factor+log_rho)
  inner_bounds = vapply(ids, function(row) mesh$inner[[row]]$bound, numeric(1L))
  inner_logtails = vapply(inner_bounds, function(bound)
    staged_tail_logs(bound, 1L)[['mass']], numeric(1L))
  log_e_mass = ls(inner_logtails)
  inner_log_tail = log_factor+log_e_mass
  inner_tail = exp(inner_log_tail)
  inner_first = setNames(numeric(p), names(theta))
  inner_first[s_names] = exp(log_factor+log(dnorm(0)*colSums(abs(prepared$w)))+log_e_mass)
  inner_first['sd_re_selection'] = exp(log_factor+log(dnorm(0)*nrow(prepared$w))+.5*log(2/pi)+log_e_mass)
  inner_first[y_names] = exp(log_factor+log(y_coefficient)+log_e_mass)
  inner_first['sd_re_outcome'] = exp(log_factor+log(sum(local))+.5*log(2/pi)+log_e_mass)
  gamma_log = vapply(seq_along(selected), function(j) log(local[j])+ls(c(
    log(2)+dnorm(inner_bounds[j], log = TRUE),
    .5*log(2/pi)+ls(inner_logtails[-j]))), numeric(1L))
  inner_first['sd_error'] = exp(log_factor+ls(gamma_log))
  inner_first['cor_re'] = exp(log_factor+log_e_mass+
    ls(c(log(abs(rho))-log(q), 2*log(1+abs(rho))-2*log(q))))
  list(nested_error = probability, nested_score_error = nested_score,
    inner_companion_log = setNames(companion, as.character(ids)), inner_outer_mass = measure,
    inner_probes = staged_inner_probes(recorder, ids),
    outer_tail = outer_tail, outer_log_tail = outer_log_tail,
    outer_tail_score = outer_first+abs(score)*outer_tail,
    inner_tail = inner_tail, inner_log_tail = inner_log_tail,
    inner_tail_score = inner_first+abs(score)*inner_tail,
    tail = outer_tail+inner_tail,
    physical_prior = list(mass = measure, first_u = prior$first_u,
      first_v = prior$first_v, second_u = prior$second_u, second_v = prior$second_v,
      absolute_uv = prior$absolute_cross,
      rho_absolute_score_bound = prior$rho_first_bound, score_moments_evaluated = TRUE),
    evidence = 'empirical finite-rule discrepancies plus analytic omitted true-model tails')
}

staged_physical_pln_panel = function(theta, prepared, control, fixed, mesh,
                                     recorder, clock, score = TRUE) {
  if (!staged_physical_pln_eligible(theta, prepared, control, fixed))
    staged_mesh_required('Retained physical PLN route no longer matches its restrictions', theta)
  options = staged_options(control)
  staged_resource_check(mesh, options, clock, staged_options_name(control))
  selected = which(prepared$d == 1)
  ids = (prepared$.staged_rows %||% seq_len(prepared$n_rows))[selected]
  p = if (score) length(theta) else 0L; rho = theta[['cor_re']]; q = 1-rho^2
  if (!length(selected) || !is.finite(q) || q <= 0)
    staged_mesh_required('Physical PLN requires observed outcomes and interior rho', theta)
  s_names = paste0('selection.', colnames(prepared$w))
  y_names = paste0('outcome.', colnames(prepared$x))
  a = drop(prepared$w %*% theta[s_names])+prepared$offset_selection
  b = drop(prepared$x %*% theta[y_names])+prepared$offset_outcome
  rule = list(quadrature = control$observation_quadrature, meshes = control$integration_mesh$inner,
    row_ids = ids, options = staged_options(control, TRUE), reference = FALSE, paired = TRUE,
    recorder = recorder, theta = theta, clock = clock, budget = control$.staged_inner_budget,
    # In physical coordinates B=b+lambda*V. Rho affects only the prior density.
    projection = list(a = NULL, b = b[selected], delta = 0,
      lambda = theta[['sd_re_outcome']], lambda_name = 'sd_re_outcome', rho = 0,
      alpha = integer(), beta = match(y_names, names(theta)),
      w = matrix(numeric(), length(selected), 0L),
      x = prepared$x[selected, , drop = FALSE]))
  axes = function(reference) list(u = staged_cell_rule(mesh$u, reference = reference),
    v = staged_cell_rule(mesh$v, reference = reference))
  canonical_axes = axes(FALSE); reference_axes = axes(TRUE)
  s_cache = y_cache = list(); outer_work = 0
  selection = function(axis, key) {
    if (!is.null(s_cache[[key]])) return(s_cache[[key]])
    U = axis$nodes[, 1L]
    value = split_probit_values(a, 2*prepared$d-1, theta[['sd_re_selection']], U, score)
    gradient = matrix(0, length(U), p)
    if (score) {
      gradient[, match(s_names, names(theta))] = value$partials %*% prepared$w
      gradient[, match('sd_re_selection', names(theta))] = U*value$sum_partials
    }
    ans = list(value = value$log_product, gradient = gradient)
    s_cache[[key]] <<- ans; ans
  }
  outcome = function(axis, key, reference) {
    canonical_key = paste(key, FALSE, sep = ':'); reference_key = paste(key, TRUE, sep = ':')
    requested_key = if (reference) reference_key else canonical_key
    if (!is.null(y_cache[[requested_key]])) return(y_cache[[requested_key]])
    V = axis$nodes[, 1L]; n = length(V); count = length(selected)
    row = rep(selected, n); input = rep(b[selected], n)+
      theta[['sd_re_outcome']]*rep(V, each = count)
    term = integrate_pln_inner_staged(NULL, input, prepared$y[row], theta[['sd_error']], 0,
      as.integer(score), lgamma(prepared$y[row]+1), rule, pln_re_inner_kernel, TRUE)
    for (key in c(canonical_key, reference_key)) {
      block = if (key == canonical_key) term else term$reference_block
      value = colSums(matrix(block$value, count, n)); gradient = matrix(0, n, p)
      if (score) {
        partial = matrix(block$gradient[, 1L], count, n)
        gradient[, match(y_names, names(theta))] = t(partial) %*% prepared$x[selected, , drop = FALSE]
        gradient[, match('sd_re_outcome', names(theta))] = V*colSums(partial)
        gradient[, match('sd_error', names(theta))] = colSums(matrix(block$gradient[, 2L], count, n))
      }
      y_cache[[key]] <<- list(value = value, gradient = gradient)
    }
    y_cache[[requested_key]]
  }
  combine = function(u_ref, v_ref, inner_ref = FALSE) {
    axis = list(u = if (u_ref) reference_axes$u else canonical_axes$u,
      v = if (v_ref) reference_axes$v else canonical_axes$v)
    u_key = if (u_ref) 'reference' else 'canonical'
    v_key = if (v_ref) 'reference' else 'canonical'
    s = selection(axis$u, u_key); y = outcome(axis$v, v_key, inner_ref)
    U = axis$u$nodes[, 1L]; V = axis$v$nodes[, 1L]
    outer_work <<- outer_work+length(U)*length(V)
    sa = s$gradient; sb = y$gradient; cross_score = numeric(p)
    constant = numeric(p)
    if (score) {
      ir = match('cor_re', names(theta))
      sa[, ir] = sa[, ir]-rho*U^2/q^2; sb[, ir] = sb[, ir]-rho*V^2/q^2
      cross_score[ir] = (1+rho^2)/q^2; constant[ir] = rho/q
    }
    cells = bilinear_cells_cpp(axis$u$log_weights-U^2/(2*q)+s$value,
      axis$v$log_weights-V^2/(2*q)+y$value, U, V, rho/q, sa, sb, cross_score,
      axis$u$cells, axis$v$cells, nrow(mesh$u$lower), nrow(mesh$v$lower))
    cells[, 1L] = cells[, 1L]-log(2*pi)-.5*log(q)
    if (score) {
      cells[, -1L] = sweep(cells[, -1L, drop = FALSE], 2L, constant, '+')
      if (length(fixed)) cells[, 1L+match(names(fixed), names(theta))] = 0
    }
    total = logsum_score(cells[, 1L], cells[, -1L, drop = FALSE], 1L)[1L, ]
    if (any(!is.finite(total)))
      staged_mesh_required('Physical PLN likelihood or finite score is not finite', theta)
    clock()
    list(cells = cells, value = total[1L], score = if (score) total[-1L] else numeric(length(theta)))
  }
  canonical = combine(FALSE, FALSE)
  comparisons = list(outer = combine(TRUE, TRUE), inner = combine(FALSE, FALSE, TRUE),
    complete = combine(TRUE, TRUE, TRUE), u = combine(TRUE, FALSE), v = combine(FALSE, TRUE))
  nc = nrow(canonical$cells); weight = exp(canonical$cells[, 1L]-canonical$value)
  cell_error = cell_refinement = numeric(nc); dimension_error = matrix(0, nc, 2L)
  score_error = numeric(length(theta)); loglik_error = 0
  for (name in names(comparisons)) {
    comparison = comparisons[[name]]; ref_weight = exp(comparison$cells[, 1L]-canonical$value)
    probability = abs(weight-ref_weight)
    cell_error = pmax(cell_error, probability)
    loglik_error = max(loglik_error, abs(canonical$value-comparison$value))
    first = if (score) abs(canonical$cells[, -1L, drop = FALSE]*weight-
      comparison$cells[, -1L, drop = FALSE]*ref_weight) else matrix(0, nc, 0L)
    refinement = probability/options$panel_loglik_tol+if (score)
      rowSums(sweep(first, 2L, options$panel_score_tol*(staged_score_floor(options)+abs(canonical$score)), '/')) else 0
    cell_refinement = pmax(cell_refinement, refinement)
    if (name %in% c('u','v')) dimension_error[, match(name, c('u','v'))] = refinement
    if (score) score_error = pmax(score_error, staged_score_discrepancy(
      canonical$cells[, -1L, drop = FALSE], comparison$cells[, -1L, drop = FALSE],
      weight, ref_weight, canonical$score, comparison$score))
  }
  prior = staged_physical_prior(canonical_axes$u, canonical_axes$v, rho, score = score, absolute_cross = score)
  clock()
  diagnostic_mesh = mesh; diagnostic_mesh$inner = control$integration_mesh$inner
  diagnostic = staged_physical_pln_diagnostics(theta, prepared, fixed, diagnostic_mesh,
    recorder, prior, canonical$value, canonical$score, derivatives = score)
  panel_ids = prepared$.staged_panel_ids %||% 1L
  recorder$outer_mass[panel_ids] = exp(prior$value)
  recorder$outer_moment[panel_ids, 1:2] = if (score)
    exp(prior$value)*c(prior$first_u, prior$first_v) else NA_real_
  prior_work = nrow(canonical_axes$u$nodes)*nrow(canonical_axes$v$nodes)+if (score)
    sum(canonical_axes$u$nodes[, 1L] != 0)*sum(canonical_axes$v$nodes[, 1L] != 0) else 0
  list(loglik = canonical$value, score = canonical$score, reference_score = comparisons$complete$score,
    loglik_error = max(loglik_error, sum(cell_error))+diagnostic$nested_error,
    score_error = score_error+diagnostic$nested_score_error+
      diagnostic$outer_tail_score+diagnostic$inner_tail_score,
    tail = diagnostic$tail, outer_tail = diagnostic$outer_tail, inner_tail = diagnostic$inner_tail,
    outer_tail_score = diagnostic$outer_tail_score, inner_tail_score = diagnostic$inner_tail_score,
    cell_error = cell_error, cell_refinement_error = cell_refinement, dimension_error = dimension_error,
    nested_error = diagnostic$nested_error, nested_score_error = diagnostic$nested_score_error,
    inner_companion_log = diagnostic$inner_companion_log, inner_outer_mass = diagnostic$inner_outer_mass,
    inner_probes = diagnostic$inner_probes, physical_prior = diagnostic$physical_prior,
    prior_normalization_defect = abs(expm1(prior$value)),
    nodes = staged_physical_node_count(mesh), intervals = nc, dimension = 2L,
    axis_intervals = c(u = nrow(mesh$u$lower), v = nrow(mesh$v$lower)),
    canonical_axis_nodes = c(u = nrow(canonical_axes$u$nodes), v = nrow(canonical_axes$v$nodes)),
    reference_axis_nodes = c(u = nrow(reference_axes$u$nodes), v = nrow(reference_axes$v$nodes)),
    outer_rule_node_work = outer_work,
    prior_rule_node_work = prior_work, reduction_log_term_evaluations = 2*(outer_work+prior_work),
    reduction_work_scope = 'six likelihood comparison grids plus the performed prior/moment grids; two streaming passes per grid',
    score_checks_available = score,
    outer_radius = staged_outer_radius(mesh), batch = list(reduction = 'streamed Cartesian cell sums'),
    coordinate_system = 'physical standard-normal U,V with correlated prior density',
    route = 'physical selected PLN with explicit zero occasion correlation')
}


