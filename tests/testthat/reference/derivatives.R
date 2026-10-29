# Retained derivative engine for reference verification and mathematical fixtures.
# Production likelihoods use the model owners and R/model_assembly.R.
# Coordinates 1:q are latent; the remaining coordinates are natural parameters.
# Only latent-latent-parameter third derivatives are stored, never a d^3 tensor.

jet = function(value, d, q = 0L, order = 1L) {
  stopifnot(d >= 0L, d == as.integer(d), q >= 0L, q <= d,
            q == as.integer(q), order %in% 0:3)
  d = as.integer(d); q = as.integer(q); order = as.integer(order)
  n = length(value)
  structure(list(value = as.numeric(value), gradient = matrix(0, n, d),
    hessian = if (order >= 2L) array(0, c(n, q, d)) else NULL,
    third = if (order >= 3L) array(0, c(n, q, q, d)) else NULL,
    d = d, q = q, order = order), class = 'ps_jet')
}

jet_subset = function(x, rows) {
  out = x
  out$value = x$value[rows]
  for (field in c('gradient', 'hessian', 'third')) {
    if (is.null(x[[field]])) next
    dimensions = dim(x[[field]])
    values = matrix(x[[field]], length(x$value), prod(dimensions[-1L]))
    out[[field]] = array(values[rows, , drop = FALSE], c(length(out$value), dimensions[-1L]))
  }
  out
}

jet_promote = function(x, context, n) {
  if (!inherits(x, 'ps_jet')) {
    if (!is.numeric(x) || !length(x) %in% c(1L, n))
      stop('Numeric derivative inputs must have matching rows or one row')
    return(jet(rep_len(x, n), context$d, context$q, context$order))
  }
  if (!identical(c(x$d, x$q, x$order), c(context$d, context$q, context$order)))
    stop('Derivative coordinates and orders must agree')
  if (length(x$value) == n) return(x)
  if (length(x$value) == 1L) return(jet_subset(x, rep(1L, n)))
  stop('Derivative inputs must have matching rows or one row')
}

jet_pair = function(x, y) {
  context = if (inherits(x, 'ps_jet')) x else y
  size = function(z) if (inherits(z, 'ps_jet')) length(z$value) else length(z)
  n = max(size(x), size(y))
  list(jet_promote(x, context, n), jet_promote(y, context, n))
}

jet_replace = function(x, rows, value) {
  indices = seq_along(x$value)[rows]
  replacement = jet_promote(value, x, length(indices))
  x$value[indices] = replacement$value
  for (field in c('gradient', 'hessian', 'third')) {
    if (is.null(x[[field]])) next
    dimensions = dim(x[[field]])
    values = matrix(x[[field]], length(x$value), prod(dimensions[-1L]))
    values[indices, ] = matrix(replacement[[field]], length(indices), ncol(values))
    x[[field]] = array(values, dimensions)
  }
  x
}

jet_where = function(condition, yes, no) {
  pair = jet_pair(yes, no)
  n = max(length(condition), length(pair[[1]]$value))
  pair = lapply(pair, jet_promote, context = pair[[1]], n = n)
  if (length(condition) == 1L) condition = rep(condition, length(pair[[1]]$value))
  if (length(condition) != length(pair[[1]]$value) || anyNA(condition))
    stop('Invalid derivative selection condition')
  jet_replace(pair[[2]], which(condition), jet_subset(pair[[1]], which(condition)))
}

jet_hessian_slice = function(x, j) {
  matrix(x$hessian[, j, , drop = FALSE], length(x$value), x$d)
}

jet_third_slice = function(x, j, k) {
  matrix(x$third[, j, k, , drop = FALSE], length(x$value), x$d)
}

jet_unary = function(x, value, first, second = 0, third = 0) {
  out = jet(value, x$d, x$q, x$order)
  if (x$order == 0L) return(out)
  out$gradient = x$gradient * first
  if (x$order >= 2L) for (j in seq_len(x$q)) {
    out$hessian[, j, ] = jet_hessian_slice(x, j)*first +
      x$gradient * (second*x$gradient[, j])
  }
  if (x$order >= 3L) for (j in seq_len(x$q)) for (k in seq_len(x$q)) {
    hj = jet_hessian_slice(x, j); hk = jet_hessian_slice(x, k)
    out$third[, j, k, ] = jet_third_slice(x, j, k)*first +
      (x$gradient*hj[, k] + hj*x$gradient[, k] + hk*x$gradient[, j])*second +
      x$gradient*(third*x$gradient[, j]*x$gradient[, k])
  }
  out
}

jet_product = function(x, y) {
  out = jet(x$value*y$value, x$d, x$q, x$order)
  if (x$order == 0L) return(out)
  out$gradient = x$gradient*y$value + y$gradient*x$value
  if (x$order >= 2L) for (j in seq_len(x$q)) {
    out$hessian[, j, ] = jet_hessian_slice(x, j)*y$value +
      jet_hessian_slice(y, j)*x$value +
      y$gradient*x$gradient[, j] + x$gradient*y$gradient[, j]
  }
  if (x$order >= 3L) for (j in seq_len(x$q)) for (k in seq_len(x$q)) {
    xj = jet_hessian_slice(x, j); xk = jet_hessian_slice(x, k)
    yj = jet_hessian_slice(y, j); yk = jet_hessian_slice(y, k)
    out$third[, j, k, ] = jet_third_slice(x, j, k)*y$value +
      jet_third_slice(y, j, k)*x$value +
      y$gradient*xj[, k] + xj*y$gradient[, k] + xk*y$gradient[, j] +
      x$gradient*yj[, k] + yj*x$gradient[, k] + yk*x$gradient[, j]
  }
  out
}

jet_power = function(x, exponent) {
  if (length(exponent) != 1L) stop('A numeric exponent must be scalar')
  if (exponent == 0) return(jet(rep(1, length(x$value)), x$d, x$q, x$order))
  if (exponent == 1) return(x)
  v = x$value
  third = if (exponent == 2) 0 else exponent*(exponent-1)*(exponent-2)*v^(exponent-3)
  jet_unary(x, v^exponent, exponent*v^(exponent-1),
    exponent*(exponent-1)*v^(exponent-2), third)
}

#' @export
#' @noRd
Ops.ps_jet = function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == '+') return(e1)
    if (.Generic == '-') return(jet_unary(e1, -e1$value, -1))
  }
  if (.Generic == '^' && !inherits(e2, 'ps_jet')) return(jet_power(e1, e2))
  pair = jet_pair(e1, e2); x = pair[[1]]; y = pair[[2]]
  if (.Generic %in% c('+', '-')) {
    sign = if (.Generic == '+') 1 else -1
    out = x
    for (field in c('value', 'gradient', 'hessian', 'third'))
      if (!is.null(out[[field]])) out[[field]] = x[[field]] + sign*y[[field]]
    return(out)
  }
  switch(.Generic, '*' = jet_product(x, y), '/' = jet_product(x, jet_power(y, -1)),
    '^' = exp(y*log(x)), stop('Unsupported derivative operator: ', .Generic))
}

#' @export
#' @noRd
Math.ps_jet = function(x, ...) {
  v = x$value
  if (.Generic == 'exp') return(jet_unary(x, exp(v), exp(v), exp(v), exp(v)))
  if (.Generic == 'log') return(jet_unary(x, log(v), 1/v, -1/v^2, 2/v^3))
  if (.Generic == 'sqrt') return(jet_power(x, 0.5))
  if (.Generic == 'tanh') {
    t = tanh(v); first = 1-t*t
    return(jet_unary(x, t, first, -2*t*first, -2*first*(1-3*t*t)))
  }
  stop('Unsupported derivative function: ', .Generic)
}

jet_log_pnorm = function(x) {
  value = pnorm(x$value, log.p = TRUE)
  mills = exp(dnorm(x$value, log = TRUE)-value)
  second = -mills*(x$value+mills)
  third = mills*((x$value+mills)*(x$value+2*mills)-1)
  tail = !is.na(x$value) & x$value < -30
  if (any(tail)) {
    t = -x$value[tail]
    mills[tail] = t+1/t-2/t^3+10/t^5-74/t^7+706/t^9-8162/t^11
    second[tail] = -1+1/t^2-6/t^4+50/t^6-518/t^8+6354/t^10-89782/t^12
    third[tail] = 2/t^3-24/t^5+300/t^7-4144/t^9+63540/t^11-1077384/t^13
  }
  jet_unary(x, value, mills, second, third)
}

jet_order = function(x, order) {
  x$order = as.integer(order)
  if (order < 2L) x$hessian = NULL
  if (order < 3L) x$third = NULL
  x
}

# Fa/F, Fb/F and Fr/F are scaled before differentiation, avoiding 0/0 tails.
# Holding log(F) fixed here gives derivatives of each CDF partial divided by F.
jet_bvn_ratios = function(a, b, r, log_value, order) {
  a = jet_order(a, order); b = jet_order(b, order); r = jet_order(r, order)
  variance = 1-r*r
  scale = sqrt(variance)
  fa = exp(-a*a/2-log(2*pi)/2 + jet_log_pnorm((b-r*a)/scale)-log_value)
  fb = exp(-b*b/2-log(2*pi)/2 + jet_log_pnorm((a-r*b)/scale)-log_value)
  fr = exp(-log(2*pi)-log(variance)/2-(a*a-2*r*a*b+b*b)/(2*variance)-log_value)
  list(fa, fb, fr)
}

jet_log_bvn = function(a, b, r) {
  inputs = list(a, b, r)
  is_jet = vapply(inputs, inherits, logical(1), 'ps_jet')
  if (!any(is_jet)) stop('A bivariate derivative input is required')
  context = inputs[[which(is_jet)[1L]]]
  n = max(vapply(inputs, function(x) if (inherits(x, 'ps_jet')) length(x$value) else length(x), integer(1)))
  inputs = lapply(inputs, jet_promote, context = context, n = n)
  a = inputs[[1L]]; b = inputs[[2L]]; r = inputs[[3L]]
  value = log_bvn_value(a$value, b$value, r$value)
  out = jet(value, a$d, a$q, a$order)
  if (a$order == 0L) return(out)
  ratios = jet_bvn_ratios(a, b, r, value, a$order-1L)
  for (i in 1:3) out$gradient = out$gradient + inputs[[i]]$gradient*ratios[[i]]$value
  if (a$order >= 2L) for (j in seq_len(a$q)) {
    raw = matrix(0, n, a$d)
    for (i in 1:3) raw = raw + jet_hessian_slice(inputs[[i]], j)*ratios[[i]]$value +
      inputs[[i]]$gradient*ratios[[i]]$gradient[, j]
    out$hessian[, j, ] = raw-out$gradient*out$gradient[, j]
  }
  if (a$order >= 3L) for (j in seq_len(a$q)) for (k in seq_len(a$q)) {
    raw = matrix(0, n, a$d)
    for (i in 1:3) {
      x = inputs[[i]]; ratio = ratios[[i]]
      raw = raw + jet_third_slice(x, j, k)*ratio$value +
        jet_hessian_slice(x, j)*ratio$gradient[, k] +
        jet_hessian_slice(x, k)*ratio$gradient[, j] +
        x$gradient*jet_hessian_slice(ratio, j)[, k]
    }
    hj = jet_hessian_slice(out, j); hk = jet_hessian_slice(out, k)
    out$third[, j, k, ] = raw - out$gradient*hj[, k] -
      hj*out$gradient[, k] - hk*out$gradient[, j] -
      out$gradient*(out$gradient[, j]*out$gradient[, k])
  }
  out
}

jet_logadd = function(x, y) {
  pair = jet_pair(x, y); x = pair[[1]]; y = pair[[2]]
  maximum = pmax(x$value, y$value)
  out = jet(maximum, x$d, x$q, x$order)
  # An absent component contributes no derivative, even when its log score is
  # infinite. Copy the surviving component before doing any 0 * Inf arithmetic.
  left = which(is.finite(x$value) & y$value == -Inf)
  right = which(is.finite(y$value) & x$value == -Inf)
  out = jet_replace(out, left, jet_subset(x, left))
  out = jet_replace(out, right, jet_subset(y, right))
  rows = which(is.finite(x$value) & is.finite(y$value))
  if (!length(rows)) return(out)
  shift = maximum[rows]
  value = log(exp(jet_subset(x, rows)-shift)+exp(jet_subset(y, rows)-shift))+shift
  jet_replace(out, rows, value)
}

jet_group_sum = function(x, group, n_groups) {
  if (length(group) != length(x$value) || anyNA(group) ||
      any(group < 1L | group > n_groups | group != as.integer(group)))
    stop('Invalid derivative group indices')
  grouped = function(value) {
    value = matrix(value, length(x$value))
    out = matrix(0, n_groups, ncol(value))
    if (!nrow(value) || !ncol(value)) return(out)
    totals = rowsum(value, group, reorder = FALSE)
    out[as.integer(rownames(totals)), ] = totals
    out
  }
  out = jet(drop(grouped(x$value)), x$d, x$q, x$order)
  for (field in c('gradient', 'hessian', 'third')) {
    if (is.null(x[[field]])) next
    out[[field]] = array(grouped(x[[field]]), c(n_groups, dim(x[[field]])[-1L]))
  }
  out
}
