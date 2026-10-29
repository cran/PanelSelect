# Independent base-R references for small, moderate-parameter mathematical checks.
# These deliberately use no production likelihood, gradient, or quadrature code.

ref_normal_rule = function(n = 15L) {
  jacobi = matrix(0, n, n)
  if (n > 1L) jacobi[cbind(seq_len(n - 1L), 2:n)] = sqrt(seq_len(n - 1L))
  eig = eigen(jacobi + t(jacobi), symmetric = TRUE)
  list(nodes = eig$values, weights = eig$vectors[1L, ]^2)
}

ref_normal_integral = function(fun, tolerance = 1e-10) {
  integrate(function(z) vapply(z, fun, numeric(1)) * dnorm(z),
            -Inf, Inf, rel.tol = tolerance, subdivisions = 300L)$value
}

# Plackett's correlation integral is independent of the production CDF routine.
ref_bvn = function(a, b, correlation) {
  base = pnorm(a) * pnorm(b)
  if (correlation == 0) return(base)
  correction = integrate(function(r) {
    exp(-(a*a - 2*r*a*b + b*b) / (2*(1-r*r))) / sqrt(1-r*r)
  }, 0, correlation, rel.tol = 1e-11)$value / (2*pi)
  base + correction
}

ref_average = function(values, scores, weights) {
  mass = values * weights
  value = sum(mass)
  list(value = value, score = drop(crossprod(mass / value, scores)))
}

# Conditional observation mass/density and its natural local score.
# Local score coordinates are selection index a, outcome index b, sd_error, tau.
ref_observation = function(family, y, d = NULL, a = 0, b = 0,
                            sd_error = 0, tau = 0, inner = ref_normal_rule()) {
  score = c(a = 0, b = 0, sd_error = 0, tau = 0)
  if (!is.null(d) && d == 0) {
    value = pnorm(-a)
    score['a'] = -dnorm(a) / value
    return(list(value = value, score = score))
  }
  selection = !is.null(d)
  root = sqrt(1 - tau*tau)
  if (family == 'probit') {
    q = 2*y - 1
    if (!selection) {
      value = pnorm(q*b)
      score['b'] = q * dnorm(b) / value
    } else {
      value = ref_bvn(a, q*b, q*tau)
      score['a'] = dnorm(a) * pnorm(q*(b-tau*a)/root) / value
      score['b'] = q * dnorm(b) * pnorm((a-tau*b)/root) / value
      density = exp(-(a*a-2*tau*a*b+b*b)/(2*root*root)) / (2*pi*root)
      score['tau'] = q * density / value
    }
  } else if (family == 'gaussian') {
    e = (y-b) / sd_error
    value = dnorm(e) / sd_error
    score['b'] = e / sd_error
    score['sd_error'] = (e*e-1) / sd_error
    if (selection) {
      index = (a+tau*e) / root
      mills = dnorm(index) / pnorm(index)
      value = value * pnorm(index)
      score['a'] = mills / root
      score['b'] = score['b'] - tau*mills/(root*sd_error)
      score['sd_error'] = score['sd_error'] - tau*e*mills/(root*sd_error)
      score['tau'] = mills * (e+tau*a) / root^3
    }
  } else if (family == 'poisson') {
    mean = exp(b)
    value = dpois(y, mean)
    score['b'] = y - mean
    if (selection) {
      value = value * pnorm(a)
      score['a'] = dnorm(a) / pnorm(a)
    }
  } else if (family == 'pln') {
    e = inner$nodes
    mean = exp(b + sd_error*e)
    mass = dpois(y, mean)
    scores = cbind(a = 0, b = y-mean, sd_error = e*(y-mean), tau = 0)
    if (selection) {
      index = (a+tau*e) / root
      mills = dnorm(index) / pnorm(index)
      mass = mass * pnorm(index)
      scores[, 'a'] = mills / root
      scores[, 'tau'] = mills * (e+tau*a) / root^3
    }
    return(ref_average(mass, scores, inner$weights))
  } else stop('Unknown reference family')
  list(value = value, score = score)
}

ref_model = function(model) {
  selected = startsWith(model, 'selection_')
  family = sub('^(selection_|.*_re$)', '', model)
  if (!selected) family = sub('_re$', '', model)
  list(family = family, selection = selected)
}

ref_parameter = function(par, name, default = 0) {
  if (name %in% names(par)) unname(par[name]) else default
}

# A panel has x, w design matrices, y, and (for selection models) d.
# Random coordinates are u and v, independently standard normal.
ref_conditional_panel = function(model, par, panel, u, v = 0,
                                  inner = ref_normal_rule()) {
  spec = ref_model(model)
  beta = par[grep('^beta', names(par))]
  alpha = par[grep('^alpha', names(par))]
  delta = ref_parameter(par, 'sd_re_selection')
  sigma = ref_parameter(par, 'sd_re_outcome')
  error = ref_parameter(par, 'sd_error')
  rho = ref_parameter(par, 'cor_re')
  tau = ref_parameter(par, 'cor_error')
  h = if (spec$selection) rho*u + sqrt(1-rho*rho)*v else u
  b = drop(panel$x %*% beta) + sigma*h
  a = if (spec$selection) drop(panel$w %*% alpha) + delta*u else rep(0, length(b))
  value = 1
  score = setNames(numeric(length(par)), names(par))
  for (j in seq_along(b)) {
    obs = ref_observation(spec$family, panel$y[j],
                           if (spec$selection) panel$d[j] else NULL,
                           a[j], b[j], error, tau, inner)
    value = value * obs$value
    score[names(beta)] = score[names(beta)] + obs$score['b'] * panel$x[j, ]
    score['sd_re_outcome'] = score['sd_re_outcome'] + obs$score['b'] * h
    if (spec$selection) {
      score[names(alpha)] = score[names(alpha)] + obs$score['a'] * panel$w[j, ]
      score['sd_re_selection'] = score['sd_re_selection'] + obs$score['a'] * u
      score['cor_re'] = score['cor_re'] + obs$score['b'] * sigma *
        (u-rho*v/sqrt(1-rho*rho))
    }
    if ('sd_error' %in% names(par)) score['sd_error'] = score['sd_error'] + obs$score['sd_error']
    if ('cor_error' %in% names(par)) score['cor_error'] = score['cor_error'] + obs$score['tau']
  }
  list(value = value, score = score[names(par)])
}

ref_panel = function(model, par, panel, outer = 15L, inner = 25L) {
  spec = ref_model(model)
  rule = ref_normal_rule(outer)
  error_rule = ref_normal_rule(inner)
  grid = if (spec$selection) expand.grid(u = seq_len(outer), v = seq_len(outer)) else
    data.frame(u = seq_len(outer), v = rep(1L, outer))
  weights = rule$weights[grid$u]
  if (spec$selection) weights = weights * rule$weights[grid$v]
  result = lapply(seq_len(nrow(grid)), function(k) {
    ref_conditional_panel(model, par, panel, rule$nodes[grid$u[k]],
                          if (spec$selection) rule$nodes[grid$v[k]] else 0,
                          error_rule)
  })
  average = ref_average(vapply(result, `[[`, numeric(1), 'value'),
                         do.call(rbind, lapply(result, `[[`, 'score')), weights)
  list(loglik = log(average$value), score = setNames(average$score, names(par)))
}

ref_difference = function(fun, par, step = 1e-5) {
  vapply(seq_along(par), function(j) {
    h = step * max(1, abs(par[j]))
    plus = minus = par
    plus[j] = plus[j] + h
    minus[j] = minus[j] - h
    (fun(plus) - fun(minus)) / (2*h)
  }, numeric(1))
}

ref_encode = function(par) {
  z = par
  scales = grepl('^sd_', names(par))
  correlations = grepl('^cor_', names(par))
  z[scales] = log(par[scales])
  z[correlations] = atanh(par[correlations])
  z
}

ref_decode = function(z) {
  par = z
  scales = grepl('^sd_', names(z))
  correlations = grepl('^cor_', names(z))
  par[scales] = exp(z[scales])
  par[correlations] = tanh(z[correlations])
  par
}

# Population mean identities derived by Gaussian integration/exponential tilting.
ref_prediction = function(model, par, x, w = NULL) {
  spec = ref_model(model)
  b = sum(x * par[grep('^beta', names(par))])
  sigma = ref_parameter(par, 'sd_re_outcome')
  error = ref_parameter(par, 'sd_error')
  response = switch(spec$family, probit = pnorm(b/sqrt(1+sigma*sigma)),
                     gaussian = b, poisson = exp(b+sigma*sigma/2),
                     pln = exp(b+(sigma*sigma+error*error)/2))
  if (!spec$selection) return(c(response = response))
  a = sum(w * par[grep('^alpha', names(par))])
  delta = ref_parameter(par, 'sd_re_selection')
  rho = ref_parameter(par, 'cor_re')
  tau = ref_parameter(par, 'cor_error')
  selection_scale = sqrt(1+delta*delta)
  selection = pnorm(a/selection_scale)
  covariance = delta*sigma*rho + if (spec$family == 'poisson') 0 else
    if (spec$family == 'probit') tau else error*tau
  joint = switch(spec$family,
    probit = ref_bvn(a/selection_scale, b/sqrt(1+sigma*sigma),
                     covariance/(selection_scale*sqrt(1+sigma*sigma))),
    gaussian = b*selection + covariance/selection_scale*dnorm(a/selection_scale),
    poisson = response*pnorm((a+covariance)/selection_scale),
    pln = response*pnorm((a+covariance)/selection_scale))
  c(response = response, selection = selection, joint = joint,
    selected = joint/selection)
}

# Integrate conditional moments directly, independently of the marginal identities.
# For a target with an independent observation error, history only reweights u/v.
ref_prediction_integral = function(model, par, x, w = NULL, history = NULL,
                                     outer = 25L, inner = 35L) {
  spec = ref_model(model)
  rule = ref_normal_rule(outer)
  error_rule = ref_normal_rule(inner)
  grid = if (spec$selection) expand.grid(u = seq_len(outer), v = seq_len(outer)) else
    data.frame(u = seq_len(outer), v = rep(1L, outer))
  weights = rule$weights[grid$u]
  if (spec$selection) weights = weights * rule$weights[grid$v]
  sigma = ref_parameter(par, 'sd_re_outcome')
  delta = ref_parameter(par, 'sd_re_selection')
  error = ref_parameter(par, 'sd_error')
  rho = ref_parameter(par, 'cor_re')
  tau = ref_parameter(par, 'cor_error')
  xb = sum(x * par[grep('^beta', names(par))])
  wa = if (spec$selection) sum(w * par[grep('^alpha', names(par))]) else 0
  moments = lapply(seq_len(nrow(grid)), function(k) {
    u = rule$nodes[grid$u[k]]
    v = if (spec$selection) rule$nodes[grid$v[k]] else 0
    b = xb + sigma * if (spec$selection) rho*u+sqrt(1-rho*rho)*v else u
    a = wa + delta*u
    response = switch(spec$family, probit = pnorm(b), gaussian = b,
                        poisson = exp(b), pln = exp(b+error*error/2))
    selection = pnorm(a)
    joint = if (!spec$selection) response else switch(spec$family,
      probit = ref_bvn(a, b, tau),
      gaussian = b*selection + error*tau*dnorm(a),
      poisson = exp(b)*selection,
      pln = sum(error_rule$weights * exp(b+error*error_rule$nodes) *
        pnorm((a+tau*error_rule$nodes)/sqrt(1-tau*tau))))
    if (!is.null(history)) weights[k] <<- weights[k] *
      ref_conditional_panel(model, par, history, u, v, error_rule)$value
    c(response = response, selection = selection, joint = joint)
  })
  average = drop(crossprod(weights/sum(weights), do.call(rbind, moments)))
  if (!spec$selection) return(average['response'])
  c(average, selected = unname(average['joint']/average['selection']))
}
