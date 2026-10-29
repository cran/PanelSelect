# Explicit PLN observation-level adaptive Gauss-Hermite integration.
# The scalar kernel is strictly log-concave. Its mode, curvature, moving nodes
# and log Jacobian are differentiated as one finite objective. Compact degree-
# three Taylor coefficients supply the mixed derivatives needed by outer AGHQ.
# This opt-in scalar R route prioritizes a complete derivative contract;
# large-panel performance and node adequacy must be assessed separately.
# Local inputs are p = (a, b, gamma, tau): the selection and outcome indices,
# the PLN observation SD (sd_error), and the wave-error correlation. The
# indices already include the current person-effect contribution.

# 1. Local Taylor algebra --------------------------------------------------
# Entries are Taylor coefficients in (a, b, gamma, tau), not raw derivatives.
# inner_taylor_derivative restores the multi-index factorials at extraction.

# The Mills ratio and its curvature M(x + M) = -F'' share normal_log_derivatives.
pln_inner_mills = function(x) {
  derivatives = normal_log_derivatives(x, 2L)
  list(ratio = derivatives$first, curvature = -derivatives$second)
}

inner_taylor_setup = local({
  grid = as.matrix(expand.grid(a = 0:3, b = 0:3, gamma = 0:3, tau = 0:3))
  grid = grid[rowSums(grid) <= 3, , drop = FALSE]
  grid = grid[order(rowSums(grid), grid[, 1], grid[, 2], grid[, 3], grid[, 4]), , drop = FALSE]
  keys = apply(grid, 1, paste, collapse = ',')
  pairs = lapply(seq_len(nrow(grid)), function(k) {
    left = which(apply(grid, 1, function(x) all(x <= grid[k, ])))
    right = match(apply(matrix(grid[k, ], length(left), 4, byrow = TRUE) -
      grid[left, , drop = FALSE], 1, paste, collapse = ','), keys)
    cbind(left, right)
  })
  # All coefficient pairs, and the 0/1 matrix summing them into each product term.
  index = do.call(rbind, pairs)
  owner = rep(seq_along(pairs), vapply(pairs, nrow, 1L))
  product = matrix(0, length(pairs), length(owner))
  product[cbind(owner, seq_along(owner))] = 1
  list(grid = grid, pairs = pairs, keys = keys, n = nrow(grid),
    left = index[, 1], right = index[, 2], product = product)
})

inner_taylor_constant = function(x) c(x, rep(0, inner_taylor_setup$n - 1L))

inner_taylor_variable = function(x, j) {
  out = inner_taylor_constant(x)
  index = rep(0, 4)
  index[j] = 1
  out[match(paste(index, collapse = ','), inner_taylor_setup$keys)] = 1
  out
}

inner_taylor_multiply = function(a, b) {
  s = inner_taylor_setup
  # One matrix product sums every pair. A nonfinite input keeps the pairwise
  # sums: a product would spread its NaN through the zero entries to every term.
  if (all(is.finite(a)) && all(is.finite(b))) drop(s$product %*% (a[s$left] * b[s$right])) else
    vapply(s$pairs, function(p) sum(a[p[, 1]] * b[p[, 2]]), numeric(1))
}

inner_taylor_power = function(a, p) {
  stopifnot(is.finite(a[1]), a[1] > 0)
  u = a
  u[1] = 0
  u = u/a[1]
  term = inner_taylor_constant(1)
  out = term
  for (k in 1:3) {
    term = inner_taylor_multiply(term, u)
    out = out + choose(p, k)*term
  }
  a[1]^p*out
}

inner_taylor_exp = function(a) {
  u = a
  u[1] = 0
  term = inner_taylor_constant(1)
  out = term
  for (k in 1:3) {
    term = inner_taylor_multiply(term, u)
    out = out + term/factorial(k)
  }
  exp(a[1])*out
}

inner_taylor_log = function(a) {
  stopifnot(is.finite(a[1]), a[1] > 0)
  u = a
  u[1] = 0
  u = u/a[1]
  term = inner_taylor_constant(1)
  out = inner_taylor_constant(log(a[1]))
  for (k in 1:3) {
    term = inner_taylor_multiply(term, u)
    out = out + (-1)^(k + 1)*term/k
  }
  out
}

inner_taylor_compose = function(a, value, derivatives) {
  u = a
  u[1] = 0
  term = inner_taylor_constant(1)
  out = inner_taylor_constant(value)
  for (k in 1:3) {
    term = inner_taylor_multiply(term, u)
    out = out + derivatives[k]*term/factorial(k)
  }
  out
}

inner_taylor_mills_derivatives = function(z) {
  # Mills ratio M and M', ..., M''''. Below -4 the recurrence below cancels
  # large terms; the continued fraction gives d^k M/dz^k = (-1)^k d^k M/dt^k.
  if (z < -4) {
    d = mills_fraction(-z, 5L)
    return(c(normal_log_derivatives(z, 1L)$first, -d[[2L]], d[[3L]], -d[[4L]], d[[5L]]))
  }
  # r'= -z*r-r^2 recursively supplies derivatives through r''''.
  coef = numeric(5)
  coef[1] = exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))
  for (n in 0:3) coef[n + 2L] = -(z*coef[n + 1L] + (if (n > 0) coef[n] else 0) +
    sum(coef[1:(n + 1L)]*rev(coef[1:(n + 1L)])))/(n + 1L)
  coef*factorial(0:4)
}

inner_taylor_mills = function(z, shift = 0L) {
  d = inner_taylor_mills_derivatives(z[1])
  inner_taylor_compose(z, d[shift + 1L], d[(shift + 2L):(shift + 4L)])
}

inner_taylor_log_phi = function(z) {
  d = inner_taylor_mills_derivatives(z[1])
  inner_taylor_compose(z, pnorm(z[1], log.p = TRUE), d[1:3])
}

# 2. Observation kernels and scalar mode ----------------------------------
# first = d h/de and information = -d^2 h/de^2 for the inner log integrand h.
inner_taylor_kernel = function(e, parameter, y, selected) {
  a = parameter[[1]]
  b = parameter[[2]]
  gamma = parameter[[3]]
  tau = parameter[[4]]
  eta = b + inner_taylor_multiply(gamma, e)
  mu = inner_taylor_exp(eta)
  value = y*eta - mu - inner_taylor_constant(lgamma(y + 1) + .5*log(2*pi)) -
    .5*inner_taylor_multiply(e, e)
  first = inner_taylor_multiply(gamma, inner_taylor_constant(y) - mu) - e
  information = inner_taylor_constant(1) +
    inner_taylor_multiply(inner_taylor_multiply(gamma, gamma), mu)
  if (selected) {
    cc = inner_taylor_power(inner_taylor_constant(1) - inner_taylor_multiply(tau, tau), -.5)
    z = inner_taylor_multiply(a + inner_taylor_multiply(tau, e), cc)
    loading = inner_taylor_multiply(tau, cc)
    value = value + inner_taylor_log_phi(z)
    first = first + inner_taylor_multiply(loading, inner_taylor_mills(z))
    information = information -
      inner_taylor_multiply(inner_taylor_multiply(loading, loading), inner_taylor_mills(z, 1L))
  }
  list(value = value, first = first, information = information)
}

inner_taylor_scalar_kernel = function(e, p, y, selected) {
  eta = p[2] + p[3]*e
  mu = exp(eta)
  value = y*eta - mu - lgamma(y + 1) + dnorm(e, log = TRUE)
  first = p[3]*(y - mu) - e
  information = 1 + p[3]^2*mu
  if (selected) {
    cc = sqrt(1 - p[4]^2)
    z = (p[1] + p[4]*e)/cc
    d = inner_taylor_mills_derivatives(z)
    value = value + pnorm(z, log.p = TRUE)
    first = first + p[4]/cc*d[1]
    information = information - (p[4]/cc)^2*d[2]
  }
  list(value = value, first = first, information = information)
}

inner_taylor_mode = function(p, y, selected) {
  score = function(e) inner_taylor_scalar_kernel(e, p, y, selected)$first
  lo = -1
  hi = 1
  while (score(lo) < 0 && abs(lo) < 8192) lo = 2*lo
  while (score(hi) > 0 && abs(hi) < 8192) hi = 2*hi
  stats::uniroot(score, c(lo, hi), tol = 1e-12)$root
}

# 3. Higher mixed derivatives of the complete adaptive rule ----------------
inner_taylor_derivative = function(jet, indices) {
  powers = tabulate(indices, nbins = 4)
  jet[match(paste(powers, collapse = ','), inner_taylor_setup$keys)]*prod(factorial(powers))
}

pln_inner_adaptive_derivatives = function(p, y, H, selected = TRUE, rule = normal_rule(H, 1L)) {
  stopifnot(length(p) == 4, abs(p[4]) < 1, p[3] >= 0)
  parameter = lapply(1:4, function(j) inner_taylor_variable(p[j], j))
  mode0 = inner_taylor_mode(p, y, selected)
  mode = inner_taylor_constant(mode0)
  # Formal implicit Newton doubles the matched polynomial order per update.
  # Four steps cover degree three and remove the scalar root's residual.
  for (k in 1:4) {
    at = inner_taylor_kernel(mode, parameter, y, selected)
    mode = mode + inner_taylor_multiply(at$first, inner_taylor_power(at$information, -1))
  }
  at = inner_taylor_kernel(mode, parameter, y, selected)
  sd = inner_taylor_power(at$information, -.5)
  gh = list(nodes = drop(rule$nodes), log_weights_1d = rule$log_weights)
  terms = lapply(seq_len(H), function(k) {
    e = mode + gh$nodes[k]*sd
    inner_taylor_kernel(e, parameter, y, selected)$value + inner_taylor_log(sd) +
      inner_taylor_constant(gh$log_weights_1d[k] + gh$nodes[k]^2/2 + .5*log(2*pi))
  })
  peak = max(vapply(terms, function(t) t[1], numeric(1)))
  total = Reduce('+', lapply(terms, function(t) inner_taylor_exp(t - inner_taylor_constant(peak))))
  value = inner_taylor_log(total) + inner_taylor_constant(peak)
  gradient = vapply(1:4, function(j) inner_taylor_derivative(value, j), numeric(1))
  hessian = matrix(0, 4, 4)
  third = array(0, c(4, 4, 4))
  for (i in 1:4) for (j in 1:4) {
    hessian[i, j] = inner_taylor_derivative(value, c(i, j))
    for (k in 1:4) third[i, j, k] = inner_taylor_derivative(value, c(i, j, k))
  }
  list(value = value[1], gradient = gradient, hessian = hessian, third = third,
    mode = mode[1], mode_jet = mode, information = at$information[1],
    mode_residual = max(abs(at$first)), jet = value)
}

# 4. Value and first-score route -------------------------------------------
# Value/first-score requests use this explicit chain rule. The Taylor route
# supplies the higher mixed derivatives required for outer adaptation.
pln_inner_adaptive_first = function(p, y, H, selected, score = TRUE, rule = normal_rule(H, 1L)) {
  gamma = p[3]
  tau = p[4]
  cc = sqrt(1 - tau*tau)
  kernel = function(e) {
    eta = p[2] + gamma*e
    mu = exp(eta)
    value = y*eta - mu - lgamma(y + 1) + dnorm(e, log = TRUE)
    first = gamma*(y - mu) - e
    information = 1 + gamma^2*mu
    if (selected) {
      z = (p[1] + tau*e)/cc
      r = pln_inner_mills(z)
      value = value + pnorm(z, log.p = TRUE)
      first = first + tau/cc*r$ratio
      information = information + (tau/cc)^2*r$curvature
    }
    list(value = value, first = first, information = information)
  }
  lo = -1
  hi = 1
  while (kernel(lo)$first < 0 && abs(lo) < 8192) lo = 2*lo
  while (kernel(hi)$first > 0 && hi < 8192) hi = 2*hi
  m = stats::uniroot(function(e) kernel(e)$first, c(lo, hi), tol = 1e-12)$root
  I = kernel(m)$information
  sd = 1/sqrt(I)
  rule = list(nodes = drop(rule$nodes), log_weights_1d = rule$log_weights)
  e = m + sd*rule$nodes
  h = kernel(e)
  terms = rule$log_weights_1d + h$value + rule$nodes^2/2 + .5*log(2*pi) + log(sd)
  value = log_sum_exp(terms)
  if (!score) return(list(value = value))

  # Differentiate the mode equation h_e = 0, then total curvature I = -h_ee.
  mu = exp(p[2] + gamma*m)
  he_j = c(0, -gamma*mu, y - mu - gamma*m*mu, 0)
  hee_j = c(0, -gamma^2*mu, -2*gamma*mu - gamma^2*m*mu, 0)
  heee = -gamma^3*mu
  if (selected) {
    z = (p[1] + tau*m)/cc
    d = inner_taylor_mills_derivatives(z)
    zt = (m + tau*p[1])/cc^3
    he_j[1] = tau/cc^2*d[2]
    he_j[4] = d[1]/cc^3 + tau/cc*d[2]*zt
    hee_j[1] = tau^2/cc^3*d[3]
    hee_j[4] = 2*tau/cc^4*d[2] + tau^2/cc^2*d[3]*zt
    heee = heee + (tau/cc)^3*d[3]
  }
  dm = he_j/I
  dI = -(hee_j + heee*dm)
  mu = exp(p[2] + gamma*e)
  hp = cbind(0, y - mu, e*(y - mu), 0)
  if (selected) {
    r = pln_inner_mills((p[1] + tau*e)/cc)$ratio
    hp[, 1] = r/cc
    hp[, 4] = r*(e + tau*p[1])/cc^3
  }
  # The H x 4 matrix combines direct scores, node motion and log-Jacobian scores.
  de = matrix(dm, H, 4, byrow = TRUE) - .5*sd*rule$nodes %o% (dI/I)
  gradient = hp + de*h$first - matrix(.5*dI/I, H, 4, byrow = TRUE)
  list(value = value, gradient = colSums(gradient*exp(terms - value)))
}

# 5. Observation rows to compact model derivatives -------------------------
# r leading coordinates are the active indices; other columns are natural
# error parameters. Only index axes are retained for higher mixed derivatives.
integrate_pln_inner_adaptive = function(A, B, y, gamma, tau, order, rule, parameters = TRUE) {
  selected = !is.null(A)
  n = length(B)
  H = length(rule$log_weights)
  coordinates = c(if (selected) 'a', 'b', if (parameters) 'error', if (selected && parameters) 'tau')
  indices = c(if (selected) 1L, 2L, if (parameters) 3L, if (selected && parameters) 4L)
  r = if (selected) 2L else 1L
  latent = if (selected) 1:2 else 2L
  out = local_block(numeric(n), coordinates, r, order)
  for (i in seq_len(n)) {
    p = c(if (selected) A[i] else 0, B[i], gamma, tau)
    a = tryCatch({
      value = if (order <= 1L) pln_inner_adaptive_first(p, y[i], H, selected, score = order == 1L, rule = rule) else
        pln_inner_adaptive_derivatives(p, y[i], H, selected, rule = rule)
      checked = c(value$value, if (order) value$gradient,
        if (order >= 2L) value$hessian, if (order >= 3L) value$third)
      if (any(!is.finite(checked)))
        integration_error('Adaptive PLN observation value or derivatives are nonfinite')
      value
    }, error = function(error) {
      integration_error(paste('Adaptive PLN observation row', i, ':', conditionMessage(error)))
    })
    out$value[i] = a$value
    if (order) out$gradient[i, ] = a$gradient[indices]
    if (order >= 2L) out$hessian[i, , ] = a$hessian[latent, indices, drop = FALSE]
    if (order >= 3L) out$third[i, , , ] = a$third[latent, latent, indices, drop = FALSE]
  }
  out
}
