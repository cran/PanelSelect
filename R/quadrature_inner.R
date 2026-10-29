# PLN observation integral.
#
# Each active PLN outcome row has its own normal observation effect e
# (only selected rows in a selection model):
#   p_t = int f_t(e) phi(e) de,  f_t(e) = Pois(y_t; exp(B_t + gamma*e)) [* Phi((A_t + tau*e)/c_tau)]
# The model file supplies the integrand as a kernel returning log f_t(e) + log w
# at a batch of (row, node) pairs, with its local partials in A, B, gamma, tau
# (pln_re_inner_kernel, probit_re_pln_re_inner_kernel). This file chooses the
# observation rule and reduces each row's nodes:
#   fixed GH (default)  integrate_pln_inner, below
#   adaptive GH         integrate_pln_inner_adaptive (quadrature_inner_adaptive.R)
#   staged GL meshes    integrate_pln_inner_staged (quadrature_staged_inner.R)
# These are row-level log probabilities and derivatives, still in local
# coordinates. Outer integration and regression-design expansion happen later;
# the two quadrature levels can use different rules.

# The observation rule of a PLN model at theta: the retained staged GL meshes
# (with the projection of the row indices that staged budgets use), or the GH
# rule, tagged 'agh' for the adaptive observation rule. point holds the model
# point's indices and scales: a, b, delta, lambda, lambda_name, rho and the
# coefficient columns alpha and beta.
observation_rule = function(theta, prepared, control, point) {
  if (identical(control$observation_quadrature, 'gl') && !isTRUE(control$.staged_seed)) {
    mesh = control$integration_mesh
    if (is.null(mesh$inner)) staged_mesh_required('PLN observation mesh is unavailable', theta)
    data = prepared$model_data
    rows = if (prepared$has_selection) data$selected$rows else data$active$rows
    ids = (prepared$.staged_rows %||% seq_len(prepared$n_rows))[rows]
    return(list(quadrature = control$observation_quadrature, meshes = mesh$inner,
      row_ids = ids, reference = isTRUE(control$.staged_reference),
      options = staged_options(control, inner = TRUE),
      recorder = control$.staged_recorder, theta = theta,
      budget = control$.staged_inner_budget,
      clock = control$.staged_clock %||% function() invisible(NULL),
      projection = list(a = point$a[rows], b = point$b[rows], delta = point$delta,
        lambda = point$lambda, lambda_name = point$lambda_name, rho = point$rho,
        alpha = point$alpha, beta = point$beta,
        w = prepared$w[rows, , drop = FALSE], x = prepared$x[rows, , drop = FALSE])))
  }
  rules = control$.rules
  inner = if (!is.null(rules) && identical(rules$observation_nodes, control$observation_nodes))
    rules$inner else normal_rule(control$observation_nodes, 1L)
  # Tag this context's copy; cached rules and old serialized GH controls
  # keep their established arithmetic. Mesh seeding needs only the maps.
  if (identical(control$observation_quadrature, 'agh')) inner$quadrature = 'agh'
  inner
}

# Fixed GH by default; adaptive GH and staged GL rules are dispatched below.
# Reduce derivatives in A, B, gamma and tau before expanding regression designs.
# The kernel belongs to its model file. It receives log GH weights so their
# addition precedes the selection term, preserving the existing arithmetic order.
integrate_pln_inner = function(A, B, y, gamma, tau, order, log_factorial, rule, kernel,
    parameters = TRUE) {
  if (identical(rule$quadrature, 'gl'))
    return(integrate_pln_inner_staged(A, B, y, gamma, tau, order,
      log_factorial, rule, kernel, parameters))
  if (identical(rule$quadrature, 'agh'))
    return(integrate_pln_inner_adaptive(A, B, y, gamma, tau, order, rule, parameters))
  selected = !is.null(A)
  coordinates = c(if (selected) 'a', 'b', if (parameters) 'error',
                  if (selected && parameters) 'tau')
  r = if (selected) 2L else 1L
  n = length(B); d = length(coordinates); nodes = length(rule$log_weights)
  out = local_block(numeric(n), coordinates, r, order)
  if (!n) return(out)
  entries = 1 + d*(order > 0L) + r*d*(order > 1L) + r*r*d*(order > 2L)
  # Half the 32 MiB target is reserved for outer maps/output. Six simultaneous
  # local blocks cover components, sums, moments and R temporaries.
  size = max(1L, floor(buffer_bytes/(8*6*nodes*entries)))
  for (first in seq.int(1L, n, by = size)) {
    rows = seq.int(first, min(n, first+size-1L))
    part = observation_mixture(kernel, A, B, y, gamma, tau, rows, rule$nodes[, 1L],
      rule$log_weights, order, log_factorial, parameters, r)
    out$value[rows] = part$value
    if (order) out$gradient[rows, ] = part$gradient
    if (order >= 2L) out$hessian[rows, , ] = part$hessian
    if (order >= 3L) out$third[rows, , , ] = part$third
  }
  out
}

# One batch of rows on one rule (nodes e, log weights): the kernel at every
# (row, node) pair, row fastest, then each row's log sum_k exp(term) and its
# derivative moments (mixture_moments; value only at order 0).
observation_mixture = function(kernel, A, B, y, gamma, tau, rows, e, log_weights, order,
                               log_factorial, parameters, r) {
  count = length(rows)
  index = rep(rows, length(log_weights))
  term = kernel(if (!is.null(A)) A[index] else NULL, B[index], y[index], gamma, tau,
    rep(e, each = count), rep(log_weights, each = count), order, log_factorial[index], parameters)
  if (!order) return(list(value = logsum_score(term$value,
    matrix(numeric(), length(index), 0L), count)[, 1L]))
  mixture_moments(term$value, term$gradient, term$hessian, term$third, n = count, q = r)
}
