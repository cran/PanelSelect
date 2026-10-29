# Scalar references for the batched staged inner recorder (production uses
# staged_record_inner_batch). They read production primitives only.

staged_inner_value = function(A, B, y, gamma, tau, order, log_factorial, mesh,
                              kernel, parameters, reference = FALSE) {
  rule = staged_inner_rule(mesh, reference)
  e = rule$nodes[, 1L]; n = length(e)
  term = kernel(if (is.null(A)) NULL else rep(A, n), rep(B, n), rep(y, n),
    gamma, tau, e, rule$log_weights, order, rep(log_factorial, n), parameters)
  r = if (is.null(A)) 1L else 2L
  if (!order) return(list(value = logsum_score(term$value,
    matrix(numeric(), n, 0L), 1L)[1L, 1L], gradient = matrix(numeric(), 1L, 0L),
    hessian = NULL, third = NULL))
  mixture_moments(term$value, term$gradient, term$hessian, term$third, n = 1L, q = r)
}

staged_inner_error_logs = function(current, previous, y, tau, mesh, selected,
                                   parameters, order) {
  logs = staged_tail_logs(mesh$bound, 1L); bounds = staged_count_bounds(y)
  probability = logsum_weights(c(staged_log_difference(current$value, previous$value),
    log(bounds$mass)+logs[['mass']]), 1L)$value
  if (!order) return(list(probability = probability, first = numeric()))
  high = max(current$value, previous$value)
  first = high+log(abs(exp(current$value-high)*drop(current$gradient)-
    exp(previous$value-high)*drop(previous$gradient)))
  c_tau = sqrt(1-tau^2)
  tails = c(if (selected) log(bounds$mass*dnorm(0)/c_tau)+logs[['mass']],
    log(bounds$first)+logs[['mass']], if (parameters) log(bounds$first)+logs[['first']],
    if (selected && parameters) log(bounds$mass)+logsum_weights(c(
      log(dnorm(0)/c_tau)+logs[['first']],
      log(abs(tau)*dnorm(1)/c_tau^2)+logs[['mass']]), 1L)$value)
  for (j in seq_along(first)) first[j] = logsum_weights(c(first[j], tails[j]), 1L)$value
  list(probability = probability, first = first)
}

staged_record_inner = function(recorder, rule, id, A, B, y, gamma, tau,
                               log_factorial, current, previous, order, parameters) {
  probe = list(row = id, A = A, B = B, y = y, log_factorial = log_factorial)
  error = staged_inner_error_logs(current, previous, y, tau, rule$meshes[[id]],
    !is.null(A), parameters, order)
  if (error$probability >= recorder$inner_error[id]) {
    recorder$inner_probes[[id]] = probe
  }
  recorder$inner_error[id] = max(recorder$inner_error[id], error$probability)
  if (!order || !length(error$first)) return(invisible(NULL))
  # PLN maps are the original independent-normal maps. Their natural first
  # pullbacks are computed here only for error accounting; the model assembler
  # still owns the actual scores and all higher mixed derivatives.
  local = match(id, rule$row_ids); base = rule$projection
  p = length(rule$theta); matrix = matrix(0, length(error$first), p)
  names_local = c(if (!is.null(A)) 'a', 'b', if (parameters) 'error',
    if (!is.null(A) && parameters) 'tau')
  b = match('b', names_local)
  if (length(base$beta)) matrix[b, base$beta] = abs(base$x[local, ])
  V = if (base$lambda > 0) (B-base$b[local])/base$lambda else 0
  lambda = match(base$lambda_name, names(rule$theta))
  if (!is.na(lambda)) matrix[b, lambda] = abs(V)
  if (!is.null(A)) {
    a = match('a', names_local)
    if (length(base$alpha)) matrix[a, base$alpha] = abs(base$w[local, ])
    U = if (base$delta > 0) (A-base$a[local])/base$delta else 0
    delta = match('sd_re_selection', names(rule$theta))
    if (!is.na(delta)) matrix[a, delta] = abs(U)
  }
  rho = match('cor_re', names(rule$theta))
  if (!is.na(rho) && base$delta > 0 && base$lambda > 0) {
    if (!is.null(A)) {
      Z = (V-base$rho*U)/sqrt(1-base$rho^2)
      matrix[b, rho] = abs(base$lambda*(U-base$rho/sqrt(1-base$rho^2)*Z))
    } else {
      # Fixed tau=0 factors A out of the inner kernel, but rho still moves B.
      # Bound its pullback by the actual outer queries (including mode search).
      matrix[b, rho] = base$lambda*recorder$outer_radius*
        (1+abs(base$rho)/sqrt(1-base$rho^2))
    }
  }
  if (parameters) {
    k = match('sd_error', names(rule$theta)); if (!is.na(k)) matrix[match('error', names_local), k] = 1
    k = match('cor_error', names(rule$theta)); if (!is.na(k) && !is.null(A))
      matrix[match('tau', names_local), k] = 1
  }
  first_error = rep(-Inf, p)
  for (k in seq_len(p)) {
    error_k = logsum_weights(error$first+log(matrix[, k]), 1L)$value
    first_error[k] = error_k
    target = match(names(rule$theta)[k], recorder$parameter_names)
    recorder$inner_score_error[id, target] = max(recorder$inner_score_error[id, target], error_k)
  }
  # A row's largest first-derivative discrepancy can occur away from its
  # largest probability discrepancy. Retain both points for shared-mesh
  # refinement, rather than repeatedly refining only the value-error point.
  candidate = max(c(-Inf, first_error[!names(rule$theta) %in% recorder$fixed_names]))
  if (candidate > -Inf && candidate >= recorder$inner_score_probe_error[id]) {
    recorder$inner_score_probes[[id]] = probe
    recorder$inner_score_probe_error[id] = candidate
  }
  invisible(NULL)
}
