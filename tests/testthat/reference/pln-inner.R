# Reference adaptive inner Gauss-Hermite quadrature and its derivatives.
# Latent derivatives and at most one structural-parameter derivative are stored.

pln_inner_normal = function(x, order, depth = 64L) {
  base = normal_log_derivatives(x, min(order, 3L))
  out = list(base$value)
  if (!order) return(out)
  out = c(out, list(base$first))
  if (order >= 2L) out[[3L]] = base$second
  if (order >= 3L) out[[4L]] = base$third
  if (order < 4L) return(out)
  m = base$first
  out[[5L]] = (-x^3+3*x)*m + (4-7*x*x)*m^2 - 12*x*m^3 - 6*m^4
  if (order >= 5L)
    out[[6L]] = (x^4-6*x*x+3)*m + (15*x^3-25*x)*m^2 +
      (50*x*x-20)*m^3 + 60*x*m^4 + 24*m^5
  tail = which(x < -4 & is.finite(x))
  if (length(tail)) {
    # Laplace's Mills continued fraction, differentiated in t = -x.
    # Reciprocal derivatives follow from u * (1/u) = 1; no large tail subtraction.
    t = -x[tail]
    u = c(list(t, rep(1, length(t))), rep(list(numeric(length(t))), order-2L))
    for (level in seq.int(depth, 1L)) {
      inverse = list(1/u[[1L]])
      for (k in seq_len(order-1L)) {
        total = numeric(length(t))
        for (j in seq_len(k)) total = total + choose(k, j)*u[[j+1L]]*inverse[[k-j+1L]]
        inverse[[k+1L]] = -inverse[[1L]]*total
      }
      u = lapply(inverse, function(value) level*value)
      u[[1L]] = u[[1L]]+t
      u[[2L]] = u[[2L]]+1
    }
    for (k in seq_len(order)) out[[k+1L]][tail] = (-1)^(k-1L)*u[[k]]
  }
  out
}

pln_inner_kernel = function(e, a, b, y, gamma, tau, log_factorial, order) {
  eta = b+gamma*e
  mean = exp(eta)
  poisson = list(y*eta-mean-log_factorial)
  if (order) poisson = c(poisson, list(y-mean), rep(list(-mean), order-1L))
  selection = !is.null(a)
  scale = sqrt(1-tau*tau)
  list(e = e, a = a, gamma = gamma, tau = tau, scale = scale,
       slope = tau/scale, poisson = poisson,
       normal = if (selection) pln_inner_normal((a+tau*e)/scale, order) else NULL)
}

# Only e derivatives, index derivatives, and at most one error/tau derivative
# occur in the compact outer-geometry blocks. This is explicit PLN calculus.
pln_inner_partial = function(kernel, e_order = 0L, coordinates = character()) {
  a = sum(coordinates == 'a'); b = sum(coordinates == 'b')
  error = sum(coordinates == 'error'); tau = sum(coordinates == 'tau')
  if (error+tau > 1L) stop('Inner partial requires at most one structural coordinate')
  n = length(kernel$e)
  value = numeric(n)
  if (!a && !tau) {
    p = kernel$poisson; gamma = kernel$gamma; degree = e_order+b
    value = gamma^e_order*p[[degree+1L]]
    if (error) {
      value = gamma^e_order*kernel$e*p[[degree+2L]]
      if (e_order) value = value + e_order*gamma^(e_order-1L)*p[[degree+1L]]
    }
  }
  if (!is.null(kernel$normal) && !b && !error) {
    f = kernel$normal; c = kernel$scale; q = kernel$slope; degree = e_order+a
    normal = q^e_order*c^(-a)*f[[degree+1L]]
    if (tau) {
      normal = q^e_order*c^(-a)*f[[degree+2L]]*(kernel$e+kernel$tau*kernel$a)/c^3
      if (e_order) normal = normal + e_order*q^(e_order-1L)*c^(-a-3L)*f[[degree+1L]]
      if (a) normal = normal + a*kernel$tau*q^e_order*c^(-a-2L)*f[[degree+1L]]
    }
    value = value+normal
  }
  if (!length(coordinates)) {
    prior = switch(as.character(e_order), '0' = -.5*kernel$e^2-.5*log(2*pi),
                    '1' = -kernel$e, '2' = rep(-1, n), numeric(n))
    value = value+prior
  }
  value
}

pln_inner_mode = function(a, b, y, gamma, tau, log_factorial, tolerance, maxit) {
  n = length(b)
  terms = function(e, rows) {
    k = pln_inner_kernel(e, if (is.null(a)) NULL else a[rows], b[rows], y[rows],
                         gamma, tau, log_factorial[rows], 2L)
    list(gradient = pln_inner_partial(k, 1L), curvature = -pln_inner_partial(k, 2L))
  }
  # A shrunk count-based start avoids enormous Poisson Newton steps at zero.
  mode = gamma*(log(y+.5)-b)/(gamma^2+1/(y+.5))
  rows = seq_len(n)
  initial = terms(mode, rows)
  if (any(!is.finite(initial$gradient)) || any(!is.finite(initial$curvature)))
    integration_error('Inner PLN starting geometry is not finite')
  lower = upper = mode
  direction = sign(initial$gradient)
  bracket = which(direction != 0)
  distance = 1
  for (iteration in seq_len(maxit)) {
    if (!length(bracket)) break
    point = mode[bracket]+direction[bracket]*distance
    gradient = terms(point, bracket)$gradient
    if (anyNA(gradient)) integration_error('Inner PLN mode bracket is invalid')
    left = direction[bracket] < 0
    lower[bracket[left]] = point[left]
    upper[bracket[!left]] = point[!left]
    bracket = bracket[gradient*direction[bracket] > 0]
    distance = 2*distance
  }
  if (length(bracket)) integration_error('Inner PLN mode could not be bracketed')
  previous_step = upper-lower
  iterations = integer(n)
  active = rows
  for (iteration in seq_len(maxit)) {
    if (!length(active)) break
    current = terms(mode[active], active)
    g = current$gradient; h = current$curvature
    if (any(!is.finite(g)) || any(!is.finite(h)) || any(h <= 0))
      integration_error('Inner PLN mode geometry is invalid')
    done = abs(g/h) <= tolerance*(1+abs(mode[active]))
    iterations[active] = iteration
    active = active[!done]; g = g[!done]; h = h[!done]
    if (!length(active)) break
    positive = g > 0
    lower[active[positive]] = mode[active[positive]]
    upper[active[!positive]] = mode[active[!positive]]
    step = g/h
    next_mode = mode[active]+step
    bisect = next_mode <= lower[active] | next_mode >= upper[active] |
      abs(step) > abs(previous_step[active])/2
    next_mode[bisect] = (lower[active[bisect]]+upper[active[bisect]])/2
    previous_step[active] = next_mode-mode[active]
    mode[active] = next_mode
  }
  if (length(active)) integration_error('Inner PLN mode did not meet its tolerance')
  final = terms(mode, rows)
  list(mode = mode, curvature = final$curvature, residual = final$gradient,
       iterations = iterations)
}

pln_inner_geometry = function(a, b, y, gamma, tau, log_factorial, order,
                             mode_tol, mode_maxit) {
  coordinates = c(if (!is.null(a)) 'a', 'b', 'error', if (!is.null(a)) 'tau')
  r = if (is.null(a)) 1L else 2L; d = length(coordinates); n = length(b)
  found = pln_inner_mode(a, b, y, gamma, tau, log_factorial, mode_tol, mode_maxit)
  kernel = pln_inner_kernel(found$mode, a, b, y, gamma, tau, log_factorial, order+2L)
  H = -pln_inner_partial(kernel, 2L)
  mode = local_block(found$mode, coordinates, r, order)
  curvature = local_block(H, coordinates, r, order)
  scale = local_block(H^(-.5), coordinates, r, order)
  jacobian = local_block(-.5*log(H), coordinates, r, order)
  if (order) for (k in seq_len(d)) {
    key = coordinates[k]
    mode$gradient[, k] = pln_inner_partial(kernel, 1L, key)/H
    curvature$gradient[, k] = -pln_inner_partial(kernel, 2L, key)-
      pln_inner_partial(kernel, 3L)*mode$gradient[, k]
    scale$gradient[, k] = -.5*H^(-1.5)*curvature$gradient[, k]
    jacobian$gradient[, k] = -.5*curvature$gradient[, k]/H
  }
  if (order >= 2L) for (j in seq_len(r)) for (k in seq_len(d)) {
    xj = coordinates[j]; xk = coordinates[k]
    mj = mode$gradient[, j]; mk = mode$gradient[, k]
    mjk = (pln_inner_partial(kernel, 1L, c(xj, xk)) +
      pln_inner_partial(kernel, 2L, xj)*mk + pln_inner_partial(kernel, 2L, xk)*mj +
      pln_inner_partial(kernel, 3L)*mj*mk)/H
    Hjk = -(pln_inner_partial(kernel, 2L, c(xj, xk)) +
      pln_inner_partial(kernel, 3L, xj)*mk + pln_inner_partial(kernel, 3L, xk)*mj +
      pln_inner_partial(kernel, 4L)*mj*mk + pln_inner_partial(kernel, 3L)*mjk)
    Hj = curvature$gradient[, j]; Hk = curvature$gradient[, k]
    mode$hessian[, j, k] = mjk
    curvature$hessian[, j, k] = Hjk
    scale$hessian[, j, k] = .75*H^(-2.5)*Hj*Hk-.5*H^(-1.5)*Hjk
    jacobian$hessian[, j, k] = -.5*(Hjk/H-Hj*Hk/H^2)
  }
  if (order >= 3L) for (j in seq_len(r)) for (l in seq_len(r)) for (k in seq_len(d)) {
    xj = coordinates[j]; xl = coordinates[l]; xk = coordinates[k]
    mj = mode$gradient[, j]; ml = mode$gradient[, l]; mk = mode$gradient[, k]
    mjl = mode$hessian[, j, l]; mjk = mode$hessian[, j, k]; mlk = mode$hessian[, l, k]
    cross = mjl*mk+mjk*ml+mlk*mj
    mjlk = (pln_inner_partial(kernel, 1L, c(xj, xl, xk)) +
      pln_inner_partial(kernel, 2L, c(xj, xl))*mk +
      pln_inner_partial(kernel, 2L, c(xj, xk))*ml +
      pln_inner_partial(kernel, 2L, c(xl, xk))*mj +
      pln_inner_partial(kernel, 2L, xj)*mlk + pln_inner_partial(kernel, 2L, xl)*mjk +
      pln_inner_partial(kernel, 2L, xk)*mjl + pln_inner_partial(kernel, 3L, xj)*ml*mk +
      pln_inner_partial(kernel, 3L, xl)*mj*mk + pln_inner_partial(kernel, 3L, xk)*mj*ml +
      pln_inner_partial(kernel, 3L)*cross + pln_inner_partial(kernel, 4L)*mj*ml*mk)/H
    Hjlk = -(pln_inner_partial(kernel, 2L, c(xj, xl, xk)) +
      pln_inner_partial(kernel, 3L, c(xj, xl))*mk +
      pln_inner_partial(kernel, 3L, c(xj, xk))*ml +
      pln_inner_partial(kernel, 3L, c(xl, xk))*mj +
      pln_inner_partial(kernel, 3L, xj)*mlk + pln_inner_partial(kernel, 3L, xl)*mjk +
      pln_inner_partial(kernel, 3L, xk)*mjl + pln_inner_partial(kernel, 4L, xj)*ml*mk +
      pln_inner_partial(kernel, 4L, xl)*mj*mk + pln_inner_partial(kernel, 4L, xk)*mj*ml +
      pln_inner_partial(kernel, 4L)*cross + pln_inner_partial(kernel, 5L)*mj*ml*mk +
      pln_inner_partial(kernel, 3L)*mjlk)
    Hj = curvature$gradient[, j]; Hl = curvature$gradient[, l]; Hk = curvature$gradient[, k]
    Hcross = curvature$hessian[, j, l]*Hk + curvature$hessian[, j, k]*Hl +
      curvature$hessian[, l, k]*Hj
    mode$third[, j, l, k] = mjlk
    curvature$third[, j, l, k] = Hjlk
    scale$third[, j, l, k] = -15/8*H^(-3.5)*Hj*Hl*Hk +
      .75*H^(-2.5)*Hcross-.5*H^(-1.5)*Hjlk
    jacobian$third[, j, l, k] = -.5*(Hjlk/H-Hcross/H^2+2*Hj*Hl*Hk/H^3)
  }
  list(mode = mode, curvature = curvature, scale = scale, jacobian = jacobian,
       diagnostics = list(max_residual = max(abs(found$residual)),
         max_scaled_residual = max(abs(found$residual/found$curvature)/(1+abs(found$mode))),
         min_curvature = min(H), max_iterations = max(found$iterations)))
}

pln_inner_node_terms = function(kernel, geometry, base, index, order) {
  coordinates = geometry$mode$coordinates
  r = geometry$mode$r; d = length(coordinates); n = length(index)
  out = local_block(pln_inner_partial(kernel), coordinates, r, order)
  if (!order) return(out)
  node = local_block(kernel$e, coordinates, r, order)
  node$gradient = geometry$mode$gradient[index, , drop = FALSE] +
    geometry$scale$gradient[index, , drop = FALSE]*base
  for (k in seq_len(d)) out$gradient[, k] = pln_inner_partial(kernel, 0L, coordinates[k]) +
    pln_inner_partial(kernel, 1L)*node$gradient[, k]
  if (order >= 2L) {
    node$hessian = geometry$mode$hessian[index, , , drop = FALSE] +
      geometry$scale$hessian[index, , , drop = FALSE]*base
    for (j in seq_len(r)) for (k in seq_len(d)) {
      xj = coordinates[j]; xk = coordinates[k]
      out$hessian[, j, k] = pln_inner_partial(kernel, 0L, c(xj, xk)) +
        pln_inner_partial(kernel, 1L, xj)*node$gradient[, k] +
        pln_inner_partial(kernel, 1L, xk)*node$gradient[, j] +
        pln_inner_partial(kernel, 2L)*node$gradient[, j]*node$gradient[, k] +
        pln_inner_partial(kernel, 1L)*node$hessian[, j, k]
    }
  }
  if (order >= 3L) {
    node$third = geometry$mode$third[index, , , , drop = FALSE] +
      geometry$scale$third[index, , , , drop = FALSE]*base
    for (j in seq_len(r)) for (l in seq_len(r)) for (k in seq_len(d)) {
      xj = coordinates[j]; xl = coordinates[l]; xk = coordinates[k]
      ej = node$gradient[, j]; el = node$gradient[, l]; ek = node$gradient[, k]
      ejl = node$hessian[, j, l]; ejk = node$hessian[, j, k]; elk = node$hessian[, l, k]
      out$third[, j, l, k] = pln_inner_partial(kernel, 0L, c(xj, xl, xk)) +
        pln_inner_partial(kernel, 1L, c(xj, xl))*ek +
        pln_inner_partial(kernel, 1L, c(xj, xk))*el +
        pln_inner_partial(kernel, 1L, c(xl, xk))*ej +
        pln_inner_partial(kernel, 1L, xj)*elk + pln_inner_partial(kernel, 1L, xl)*ejk +
        pln_inner_partial(kernel, 1L, xk)*ejl + pln_inner_partial(kernel, 2L, xj)*el*ek +
        pln_inner_partial(kernel, 2L, xl)*ej*ek + pln_inner_partial(kernel, 2L, xk)*ej*el +
        pln_inner_partial(kernel, 2L)*(ejl*ek+ejk*el+elk*ej) +
        pln_inner_partial(kernel, 3L)*ej*el*ek + pln_inner_partial(kernel, 1L)*node$third[, j, l, k]
    }
  }
  out
}

local_pln_adaptive = function(b, y, gamma, order, log_factorial, rule, a = NULL,
                             tau = 0, mode_tol = 1e-10, mode_maxit = 60L) {
  coordinates = c(if (!is.null(a)) 'a', 'b', 'error', if (!is.null(a)) 'tau')
  r = if (is.null(a)) 1L else 2L; d = length(coordinates); n = length(b)
  if (length(order) != 1L || !order %in% 0:3) stop('Inner derivative order must be 0, 1, 2 or 3')
  if (!is.finite(gamma) || gamma <= 0 || !is.finite(tau) || abs(tau) >= 1)
    integration_error('Adaptive inner PLN parameters are not interior')
  out = local_block(numeric(n), coordinates, r, order)
  if (!n) return(out)
  nodes = length(rule$log_weights)
  if (!nodes || !identical(dim(rule$nodes), c(nodes, 1L))) stop('Inner rule must be one-dimensional')
  entries = 1+d*(order > 0L)+r*d*(order > 1L)+r*r*d*(order > 2L)
  size = max(1L, floor(16*1024^2/(8*8*nodes*entries)))
  diagnostics = list(max_residual = 0, max_scaled_residual = 0, min_curvature = Inf, max_iterations = 0L)
  for (start in seq.int(1L, n, by = size)) {
    rows = seq.int(start, min(n, start+size-1L)); count = length(rows)
    selected_index = if (is.null(a)) NULL else a[rows]
    geometry = pln_inner_geometry(selected_index, b[rows], y[rows], gamma, tau,
                                 log_factorial[rows], order, mode_tol, mode_maxit)
    index = rep(seq_len(count), nodes)
    base = rep(rule$nodes[, 1L], each = count)
    e = geometry$mode$value[index]+geometry$scale$value[index]*base
    kernel = pln_inner_kernel(e, if (is.null(a)) NULL else selected_index[index],
      b[rows][index], y[rows][index], gamma, tau, log_factorial[rows][index], max(1L, order))
    term = pln_inner_node_terms(kernel, geometry, base, index, order)
    # Standard-normal weights: de = C dz and phi(z) is divided out.
    term$value = term$value + rep(rule$log_weights, each = count) + .5*base^2+.5*log(2*pi)
    reduced = mixture_moments(term$value, term$gradient, term$hessian, term$third,
                              n = count, q = if (order) r else 0L)
    out$value[rows] = reduced$value+geometry$jacobian$value
    if (order) out$gradient[rows, ] = reduced$gradient+geometry$jacobian$gradient
    if (order >= 2L) out$hessian[rows, , ] = reduced$hessian+geometry$jacobian$hessian
    if (order >= 3L) out$third[rows, , , ] = reduced$third+geometry$jacobian$third
    for (name in names(diagnostics)) diagnostics[[name]] =
      if (name == 'min_curvature') min(diagnostics[[name]], geometry$diagnostics[[name]]) else
        max(diagnostics[[name]], geometry$diagnostics[[name]])
  }
  if (any(!is.finite(out$value)) || any(!is.finite(out$gradient)) ||
      any(!is.finite(out$hessian)) || any(!is.finite(out$third)))
    integration_error('Adaptive inner PLN contribution or derivatives are not finite')
  out$inner_diagnostics = diagnostics
  out
}
