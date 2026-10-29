# Stable probability values and compact normal derivatives. Ordinary cases
# retain the existing compiled distribution routines.
normal_mills = function(x) normal_log_derivatives(x, 1L)$first

# Conditional-normal integration, centered and scaled at its constrained mode.
# Scaling is important when a rare joint event has a very narrow boundary peak.
log_bvn_tail = function(a, b, r) {
  if (a > b) { tmp = a; a = b; b = tmp }
  s = sqrt(1 - r*r)
  log_kernel = function(z) stats::dnorm(z, log = TRUE) +
    stats::pnorm((b - r*z)/s, log.p = TRUE)
  slope = function(z) -z - r/s * normal_mills((b - r*z)/s)
  mode = a
  if (slope(a) < 0) {
    width = 1
    while (slope(a - width) < 0 && width < 1e8) width = 2*width
    mode = tryCatch(stats::uniroot(slope, c(a - width, a), tol = 1e-11)$root,
      error = function(error) integration_error(paste("Bivariate normal tail mode:",
                                                     conditionMessage(error))))
  }
  v = (b - r*mode)/s
  mills = normal_mills(v)
  curvature = 1 + (r/s)^2 * mills * (v + mills)
  scale = 1/sqrt(max(1, curvature) + slope(mode)^2)
  shift = log_kernel(mode)
  log_cdf_mode = stats::pnorm(v, log.p = TRUE)
  centered = function(t) {
    delta = scale*t
    dv = -r/s*delta
    at = v + dv
    log_ratio = stats::pnorm(at, log.p = TRUE) - log_cdf_mode
    tail = at < -30 & v < -30
    # log Phi = log phi - log Mills. Compute differences before squaring;
    # subtracting two enormous negative log probabilities loses tail precision.
    if (any(tail)) log_ratio[tail] = -dv[tail]*(v + dv[tail]/2) -
      log1p((normal_mills(at[tail]) - mills)/mills)
    exp(-delta*(mode + delta/2) + log_ratio)
  }
  integral = tryCatch(stats::integrate(centered, -Inf, (a - mode)/scale,
    rel.tol = 2e-11, abs.tol = 0, subdivisions = 500L),
    error = function(error) integration_error(paste("Bivariate normal tail integration:",
                                                   conditionMessage(error))))
  if (!is.finite(integral$value) || integral$value <= 0)
    integration_error("Bivariate normal tail integration failed")
  shift + log(scale) + log(integral$value)
}

gk21_rule = function() {
  # QUADPACK dqk21, https://www.netlib.org/quadpack/dqk21.f. G10 is embedded
  # at the even positive abscissae. The two rules integrate Lebesgue measure.
  x = c(.9956571630258080807, .9739065285171717201, .9301574913557082260,
        .8650633666889845107, .7808177265864168971, .6794095682990244062,
        .5627571346686046833, .4333953941292471908, .2943928627014601981,
        .1488743389816312109)
  w = c(.0116946388673718743, .0325581623079647275, .0547558965743519960,
        .0750396748109199528, .0931254545836976055, .1093871588022976419,
        .1234919762620658511, .1347092173114733259, .1427759385770600808,
        .1477391049013384914)
  g = numeric(21L)
  g[seq.int(2L, 10L, 2L)] = c(.0666713443086881376, .1494513491505805931,
    .2190863625159820440, .2692667193099963551, .2955242247147528702)
  g[seq.int(12L, 20L, 2L)] = rev(g[seq.int(2L, 10L, 2L)])
  list(nodes = c(-x, 0, rev(x)), weights = c(w, .1494455540029169057, rev(w)),
       reference_weights = g)
}

# log_bvn_tail for many elements at once, with the same centred, scaled
# integrand g(t) (g(0) = 1). The constrained mode solves slope = 0 by safeguarded
# Newton inside [a + slope(a), a]: the slope's derivative is at most -1, so that
# bracket holds the unique root. Gauss-Kronrod 21/10 panels march outward from
# t = 0, graded and then geometric. An element stops once log-concavity bounds
# the mass beyond the current edge T by g(T)|T|/|log g(T)| <= 1e-14 of its
# integral. Elements failing the scalar rel.tol, the tail bound or the mode
# search return NA for the scalar path, as do extreme arguments whose failures
# the scalar routine reports.
log_bvn_tail_conditional = function(a, b, r, chunk = 4096L) {
  n = length(a)
  if (n > chunk) return(unlist(lapply(split(seq_len(n), ceiling(seq_len(n)/chunk)),
    function(i) log_bvn_tail_conditional(a[i], b[i], r[i], chunk))))
  low = pmin(a, b); b = pmax(a, b); a = low
  s = sqrt(1 - r*r); k = r/s
  slope = function(z, i) -z - k[i]*normal_mills((b[i] - r[i]*z)/s[i])
  mode = a; converged = rep(TRUE, n)
  inside = which(slope(a, seq_len(n)) < 0)
  if (length(inside)) {
    z = upper = a[inside]; lower = z + slope(z, inside)
    for (iteration in seq_len(100L)) {
      d = normal_log_derivatives((b[inside] - r[inside]*z)/s[inside], 2L)
      f = -z - k[inside]*d$first
      lower = ifelse(f > 0, z, lower); upper = ifelse(f < 0, z, upper)
      proposal = z + f/(1 + k[inside]^2*(-d$second))
      outside = !is.finite(proposal) | proposal < lower | proposal > upper
      proposal[outside] = (lower[outside] + upper[outside])/2
      done = abs(proposal - z) <= 1e-12*pmax(1, abs(z))
      z = proposal
      if (all(done)) break
    }
    mode[inside] = z; converged[inside] = done
  }
  v = (b - r*mode)/s
  d = normal_log_derivatives(v, 2L)
  mills = d$first; log_cdf = d$value
  scale = 1/sqrt(pmax(1, 1 + k^2*(-d$second)) + slope(mode, seq_len(n))^2)
  upper = (a - mode)/scale
  # log g at a matrix t whose rows are the elements `rows`.
  log_g = function(rows, t) {
    delta = scale[rows]*t; dv = -k[rows]*delta; at = v[rows] + dv
    ratio = stats::pnorm(at, log.p = TRUE) - log_cdf[rows]
    deep = at < -30 & v[rows] < -30
    if (any(deep)) {
      row = rows[(which(deep) - 1L) %% length(rows) + 1L]; e = dv[deep]
      ratio[deep] = -e*(v[row] + e/2) - log1p((normal_mills(at[deep]) - mills[row])/mills[row])
    }
    -delta*(mode[rows] + delta/2) + ratio
  }
  rule = gk21_rule()
  kronrod = gauss = numeric(n)
  add = function(rows, lo, hi) {
    half = (hi - lo)/2
    g = exp(log_g(rows, (lo + hi)/2 + outer(half, rule$nodes)))
    kronrod[rows] <<- kronrod[rows] + half*drop(g %*% rule$weights)
    gauss[rows] <<- gauss[rows] + half*drop(g %*% rule$reference_weights)
  }
  # Mass beyond an edge at distance T with log g(T) = h: g(T) T/|h| (h < 0).
  settled = function(rows, h, distance) h < 0 & exp(h)*distance/(-h) <= 1e-14*kronrod[rows]
  edges = c(0, 1.5, 3, 5, 8, 12, 18, 26, 36, 48*1.5^(0:20))
  # Between the mode and the boundary a (t in [0, upper]); empty for a boundary mode.
  open = which(upper > 0)
  for (j in seq_len(length(edges) - 1L)) {
    if (!length(open)) break
    hi = pmin(edges[j + 1L], upper[open])
    add(open, rep(edges[j], length(open)), hi)
    finished = hi >= upper[open] | settled(open, log_g(open, matrix(hi)), hi)
    open = open[!finished]
  }
  unresolved = open
  # Below the mode, towards minus infinity.
  open = seq_len(n)
  for (j in seq_len(length(edges) - 1L)) {
    if (!length(open)) break
    add(open, rep(-edges[j + 1L], length(open)), rep(-edges[j], length(open)))
    finished = settled(open, log_g(open, matrix(-edges[j + 1L], length(open), 1L)), edges[j + 1L])
    open = open[!finished]
  }
  unresolved = union(unresolved, open)
  value = stats::dnorm(mode, log = TRUE) + log_cdf + log(scale) + log(kronrod)
  extreme = pmax(abs(a), abs(b)) > 1e6 | abs(r) > 1 - 1e-10
  valid = converged & !extreme & is.finite(value) & kronrod > 0 &
    abs(kronrod - gauss) <= 2e-11*kronrod
  valid[unresolved] = FALSE
  value[!valid] = NA_real_
  value
}

# For a+b <= 0, Plackett's identity starts at zero probability when r = -1:
# P = (2*pi)^-1 integral_0^acos(-r) exp(-E(theta)) dtheta, where
# E = A/sin(theta/2)^2 + B/cos(theta/2)^2, A=(a+b)^2/8, B=(a-b)^2/8.
# Its analytic mode avoids the normal CDF at every quadrature node. Integrate
# each side of that mode with the same K21/G10 target and log-concave tail bound
# as the conditional route. Failed error estimates return NA for that route.
log_bvn_tail_angular = function(a, b, r) {
  n = length(a)
  A = (a+b)^2/8; B = (a-b)^2/8
  angle = acos(-r)
  peak = 2*atan2(sqrt(abs(a+b)), sqrt(abs(a-b)))
  mode = pmin(peak, angle)
  inside = peak <= angle
  x0 = sin(mode/2)^2; c0 = cos(mode/2)^2
  aa = A/x0^2; aa[A == 0] = 0
  bb = B/c0^2
  linear = 2*(b-r*a)*(r*b-a)/(1-r*r)^2
  linear[inside] = 0
  # A=0 puts the mode at theta=0, not at a stationary point in x=sin(theta/2)^2.
  linear[A == 0] = B[A == 0]
  linear_c0 = linear*c0
  curvature = (aa*x0 + bb*c0)/2 + 1.5*(aa*c0 + bb*x0)
  scale = 1/sqrt(pmax(1, curvature) + (linear*sin(mode)/2)^2)
  shift = -(aa*x0 + bb*c0)
  shift[inside] = -.5*pmax(abs(a[inside]), abs(b[inside]))^2
  # Center E algebraically: E(x)-E(x0) = dx*E'(x0) +
  # dx^2*(A/(x*x0^2) + B/((1-x)*(1-x0)^2)). Both terms are nonnegative
  # on either integration side; subtracting enormous exponents is avoided.
  log_g = function(rows, t) {
    delta = scale[rows]*t
    dx = sin(mode[rows]+delta/2)*sin(delta/2)
    x = x0[rows]+dx; c = c0[rows]-dx
    small = x <= .01*x0[rows] | c <= .01*c0[rows]
    if (any(small)) {
      theta = mode[rows]+delta
      x[small] = sin(theta[small]/2)^2
      c[small] = cos(theta[small]/2)^2
    }
    # With linear = bb-aa, the centered exponent factors as
    # dx/c * (linear*c0 + aa*dx/x); the terms retain the same signs.
    # The removable 0/0 when A=x=0 gives the limit B*x/c.
    centered = aa[rows]*dx/x; centered[is.nan(centered)] = 0
    -dx/c * (linear_c0[rows] + centered)
  }
  rule = gk21_rule()
  kronrod = error = numeric(n)
  add = function(rows, lo, hi) {
    half = (hi-lo)/2
    g = exp(log_g(rows, (lo+hi)/2 + outer(half, rule$nodes)))
    K = half*drop(g %*% rule$weights)
    G = half*drop(g %*% rule$reference_weights)
    kronrod[rows] <<- kronrod[rows] + K
    error[rows] <<- error[rows] + abs(K-G)
  }
  settled = function(rows, h, distance) h < 0 & exp(h)*distance/(-h) <= 1e-14*kronrod[rows]
  # Split the first scaled interval on each side of the mode: resolving its
  # embedded error here avoids restarting many smooth tails conditionally.
  edges = c(0, 2, 4, 8, 16, 32, 64*1.5^(0:20))
  unresolved = integer()
  for (direction in c(1, -1)) {
    boundary = if (direction == 1) (angle-mode)/scale else mode/scale
    open = which(boundary > 0)
    for (j in seq_len(length(edges)-1L)) {
      if (!length(open)) break
      hi = pmin(edges[j+1L], boundary[open])
      if (direction == 1) add(open, rep(edges[j], length(open)), hi) else
        add(open, -hi, rep(-edges[j], length(open)))
      finished = hi >= boundary[open] |
        settled(open, log_g(open, matrix(direction*hi)), hi)
      finished[is.na(finished)] = FALSE
      open = open[!finished]
    }
    unresolved = union(unresolved, open)
  }
  value = shift + log(scale) + log(kronrod) - log(2*pi)
  valid = is.finite(value) & kronrod > 0 & error <= 2e-11*kronrod
  valid[unresolved] = FALSE
  value[!valid] = NA_real_
  value
}

# Accelerate only a tested, moderate-correlation domain. Other inputs and
# angular rules that fail their error estimate keep the conditional calculation.
log_bvn_tail_batch = function(a, b, r, chunk = 4096L) {
  n = length(a)
  if (n > chunk) return(unlist(lapply(split(seq_len(n), ceiling(seq_len(n)/chunk)),
    function(i) log_bvn_tail_batch(a[i], b[i], r[i], chunk))))
  value = rep(NA_real_, n)
  safe = pmax(abs(a), abs(b)) <= 40 & abs(r) <= .95
  positive = which(safe & r > 0)
  if (length(positive)) {
    lo = pmin(a[positive], b[positive]); hi = pmax(a[positive], b[positive])
    # For r>0 and z<=lo, Phi((hi-r*z)/s)>=Phi((hi-r*lo)/s). Thus
    # Phi(lo)*Phi((hi-r*lo)/s) <= P <= Phi(lo); this is a log-error bound,
    # not a probability floor. The analytical score numerators are unchanged.
    bound = -stats::pnorm((hi-r[positive]*lo)/sqrt(1-r[positive]^2), log.p = TRUE)
    marginal = positive[bound <= 1e-14]
    value[marginal] = stats::pnorm(pmin(a[marginal], b[marginal]), log.p = TRUE)
  }
  angular = which(safe & a+b <= 0 & is.na(value))
  if (length(angular)) value[angular] = log_bvn_tail_angular(a[angular], b[angular], r[angular])
  remaining = which(is.na(value))
  if (length(remaining)) value[remaining] =
    log_bvn_tail_conditional(a[remaining], b[remaining], r[remaining], chunk)
  value
}

log_bvn_value = function(a, b, r) {
  n = max(length(a), length(b), length(r))
  if (!n) return(numeric())
  a = rep_len(a, n); b = rep_len(b, n); r = rep_len(r, n)
  if (any(!is.finite(c(a, b, r))) || any(abs(r) >= 1))
    integration_error("Bivariate normal arguments must be finite with abs(correlation) < 1")
  ans = numeric(n)
  independent = r == 0
  ans[independent] = stats::pnorm(a[independent], log.p = TRUE) +
    stats::pnorm(b[independent], log.p = TRUE)
  ix = which(!independent)
  if (length(ix)) {
    probability = pbivnorm::pbivnorm(a[ix], b[ix], r[ix])
    ans[ix] = suppressWarnings(log(probability))
    tail = ix[!is.finite(probability) | probability < 1e-9]
    if (length(tail)) {
      ans[tail] = log_bvn_tail_batch(a[tail], b[tail], r[tail])
      for (j in tail[is.na(ans[tail])]) ans[j] = log_bvn_tail(a[j], b[j], r[j])
    }
  }
  ans
}

# Fixed tensor nodes revisit identical selection indices. Keep this cache local
# to one parameter evaluation; its numeric payload is bounded independently of H.
cache_normal_derivatives = function(limit = buffer_bytes, max_entries = 32L) {
  inputs = orders = values = list(); bytes = 0
  function(x, order) {
    n = length(x); first = if (n) x[[1L]] else NULL
    for (i in seq_along(inputs))
      if (n == length(inputs[[i]]) && identical(order, orders[[i]]) &&
          identical(first, if (n) inputs[[i]][[1L]] else NULL, num.eq = FALSE) &&
          identical(x, inputs[[i]], num.eq = FALSE))
        return(values[[i]])
    value = normal_log_derivatives(x, order)
    size = 8*(length(x) + sum(lengths(value)))
    if (length(inputs) < max_entries && bytes+size <= limit) {
      i = length(inputs)+1L
      inputs[[i]] <<- x; orders[[i]] <<- order; values[[i]] <<- value
      bytes <<- bytes+size
    }
    value
  }
}

# Derivatives in t = -x of the Mills ratio M = phi(x)/Phi(x), for t > 4, from
# Laplace's continued fraction M = t + 1/(t + 2/(t + 3/(t + ...))). Each level
# is a reciprocal whose derivatives follow from u*(1/u) = 1, so no large terms
# cancel. The depth falls as t grows; it was validated against 1024-bit
# references (tests/testthat/fixtures/normal-tail.R). Returns d^k M/dt^k for
# k = 0, ..., order - 1.
mills_fraction = function(t, order) {
  out = rep(list(numeric(length(t))), order)
  depth = ifelse(t < 8, 64L, ifelse(t < 16, 32L, 16L))
  for (levels in unique(depth)) {
    i = which(depth == levels); ti = t[i]
    u = c(list(ti, rep(1, length(ti))), rep(list(numeric(length(ti))), order))[seq_len(order)]
    for (level in seq.int(levels, 1L)) {
      inverse = list(1/u[[1L]])
      for (k in seq_len(order - 1L)) {
        total = 0
        for (j in seq_len(k)) total = total + choose(k, j)*u[[j + 1L]]*inverse[[k - j + 1L]]
        inverse[[k + 1L]] = -inverse[[1L]]*total
      }
      u = lapply(inverse, `*`, level)
      u[[1L]] = u[[1L]] + ti
      if (order > 1L) u[[2L]] = u[[2L]] + 1
    }
    for (k in seq_len(order)) out[[k]][i] = u[[k]]
  }
  out
}

normal_log_derivatives = function(x, order) {
  value = stats::pnorm(x, log.p = TRUE)
  if (!order) return(list(value = value))
  first = exp(-0.5*x*x - 0.5*log(2*pi) - value)
  second = if (order >= 2L) -first*(x + first) else NULL
  third = if (order >= 3L) first*((x + first)*(x + 2*first) - 1) else NULL
  tail = !is.na(x) & x < -30
  if (any(tail)) {
    t = -x[tail]
    first[tail] = t + 1/t - 2/t^3 + 10/t^5 - 74/t^7 + 706/t^9 - 8162/t^11
  }
  # Below -4 the closed forms cancel large terms; F'' = -M' and F''' = M'' in t.
  if (order >= 2L) {
    deep = !is.na(x) & x < -4
    if (any(deep)) {
      d = mills_fraction(-x[deep], order)
      second[deep] = -d[[2L]]
      if (order >= 3L) third[deep] = d[[3L]]
    }
  }
  list(value = value, first = first, second = second, third = third)
}

# Extend an already evaluated first score with the latent/mixed blocks.
normal_log_higher = function(out, f, derivative, mixed, r, order) {
  d = ncol(derivative); n = nrow(derivative)
  if (order >= 2L) for (i in seq_len(r)) {
    mi = matrix(mixed[, i, ], n, d)
    out$hessian[, i, ] = derivative*(f$second*derivative[, i]) + mi*f$first
    if (order >= 3L) for (j in seq_len(r))
      out$third[, i, j, ] = derivative*(f$third*derivative[, i]*derivative[, j]) +
        (mi*derivative[, j] + matrix(mixed[, j, ], n, d)*derivative[, i])*f$second
  }
  out
}

# Signed normal log probabilities, including exact unselected contributions.
local_univariate_probit = function(index, sign, name, order) {
  f = normal_log_derivatives(sign*index, order)
  out = local_block(f$value, name, 1L, order, if (order) f$first*sign)
  if (order >= 2L) out$hessian[, 1L, 1L] = f$second*sign^2
  if (order >= 3L) out$third[, 1L, 1L, 1L] = f$third*sign^3
  out
}

# Poisson chain derivatives used by compact observation-effect kernels.
poisson_log_higher = function(out, derivative, mu, r, order) {
  if (order >= 2L) for (i in seq_len(r)) {
    out$hessian[, i, ] = derivative*(-mu*derivative[, i])
    if (order >= 3L) for (j in seq_len(r))
      out$third[, i, j, ] = derivative*(-mu*derivative[, i]*derivative[, j])
  }
  out
}

# Compact higher derivatives of log Phi_2; model owners supply the visible
# first partials and the signs of their threshold/correlation coordinates.
binary_log_higher = function(out, a, b, r, first, signs, order) {
  n = length(a); variance = 1-r*r
  fa = first[, 1L]; fb = first[, 2L]; fr = first[, 3L]
  if (order < 2L) return(out)
  ar = (r*b-a)/variance; br = (r*a-b)/variance
  rr = r/variance + (a*b*(1+r*r)-r*(a*a+b*b))/variance^2
  raw = array(0, c(n, 2L, 3L))
  raw[, 1L, ] = cbind(-a*fa-r*fr, fr, fr*ar)
  raw[, 2L, ] = cbind(fr, -b*fb-r*fr, fr*br)
  for (i in 1:2) out$hessian[, i, ] =
    (matrix(raw[, i, ], n, 3L) - first*first[, i])*signs*signs[, i]
  if (order < 3L) return(out)
  third = array(0, c(n, 2L, 2L, 3L))
  third[, 1L, 1L, ] = cbind((a*a-1)*fa+a*r*fr-r*fr*ar,
                            -a*fr-r*fr*br, -a*fr*ar-fr-r*fr*rr)
  third[, 1L, 2L, ] = cbind(-a*fr-r*fr*br, -b*fr-r*fr*ar, fr*rr)
  third[, 2L, 1L, ] = third[, 1L, 2L, ]
  third[, 2L, 2L, ] = cbind(-b*fr-r*fr*ar, (b*b-1)*fb+b*r*fr-r*fr*br,
                            -b*fr*br-fr-r*fr*rr)
  for (i in 1:2) for (j in 1:2)
    out$third[, i, j, ] = (matrix(third[, i, j, ], n, 3L) -
      first*raw[, i, j] - matrix(raw[, i, ], n, 3L)*first[, j] -
      matrix(raw[, j, ], n, 3L)*first[, i] + 2*first*(first[, i]*first[, j]))*
      signs*(signs[, i]*signs[, j])
  out
}
