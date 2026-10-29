# Independent mathematical references. These verify integration and analytical
# scores separately from the low-order finite-rule derivative regression suite.
test_that('split singleton values and all natural scores match an exact arcsine law', {
  p = list(w = cbind('(Intercept)' = 1, x = .7),
    x = cbind('(Intercept)' = 1, x = -.4), selection_sign = 1,
    outcome_sign = 1, offset_selection = 0, offset_outcome = 0)
  theta = c('selection.(Intercept)' = 0, selection.x = 0,
    'outcome.(Intercept)' = 0, outcome.x = 0,
    sd_re_selection = 1.2, sd_re_outcome = .8, cor_re = 0, cor_error = 0)
  options = list(bound = 12, max_nodes = 1000000L)
  for (scale in list(c(.8, 1.1), c(2.2, 9))) for (rho in c(-.6, 0, .6)) for (sign_y in c(-1, 1)) {
    th = theta; th[c('sd_re_selection', 'sd_re_outcome', 'cor_re')] = c(scale, rho)
    p$outcome_sign = sign_y
    ds = sqrt(1 + scale[1]^2); ls = sqrt(1 + scale[2]^2)
    r = sign_y*rho*prod(scale)/(ds*ls)
    probability = .25 + asin(r)/(2*pi)
    density = 1/(2*pi*sqrt(1 - r^2)); exact = th*0
    exact[1:2] = dnorm(0)*.5/ds*c(1, .7)/probability
    exact[3:4] = sign_y*dnorm(0)*.5/ls*c(1, -.4)/probability
    exact['sd_re_selection'] = density*sign_y*rho*scale[2]/(ls*ds^3)/probability
    exact['sd_re_outcome'] = density*sign_y*rho*scale[1]/(ds*ls^3)/probability
    exact['cor_re'] = density*sign_y*prod(scale)/(ds*ls)/probability
    z = split_physical_panel(th, p, 20L, options, fixed = c(cor_error = 0))
    expect_equal(z$loglik, log(probability), tolerance = 1e-9)
    expect_lt(max(abs(z$score-exact)/(1+abs(exact))), 1e-8)
    expect_identical(z$loglik, split_physical_panel(th, p, 20L, options, FALSE,
      fixed = c(cor_error = 0))$loglik)
  }
})

test_that('independent one-dimensional integrals verify rho-zero values and scores', {
  th = c('selection.(Intercept)' = .21, selection.x = -.17,
    'outcome.(Intercept)' = -.28, outcome.x = .32,
    sd_re_selection = 1.2, sd_re_outcome = .8, cor_re = 0, cor_error = 0)
  p = list(w = cbind('(Intercept)' = 1, x = c(-.8, .35, .9)),
    x = cbind('(Intercept)' = 1, x = c(-.8, .9)),
    selection_sign = c(1,-1,1), outcome_sign = c(-1,1),
    offset_selection = c(.05,-.03,.08), offset_outcome = c(-.02,.04))
  reference = function(index, signs, scale) {
    integrand = function(z) dnorm(z)*exp(rowSums(pnorm(sweep(sweep(
      outer(z, rep(scale, length(index)), '*'), 2L, index, '+'),
      2L, signs, '*'), log.p = TRUE)))
    L = integrate(integrand, -Inf, Inf, rel.tol = 1e-11, abs.tol = 1e-12)
    m = integrate(function(z) z*integrand(z), -Inf, Inf,
      rel.tol = 1e-11, abs.tol = 1e-12)
    expect_identical(L$message, 'OK'); expect_identical(m$message, 'OK')
    c(value = L$value, mean = m$value/L$value)
  }
  s = reference(drop(p$w %*% th[1:2])+p$offset_selection,
    p$selection_sign, th['sd_re_selection'])
  y = reference(drop(p$x %*% th[3:4])+p$offset_outcome,
    p$outcome_sign, th['sd_re_outcome'])
  z = split_physical_panel(th, p, 20L, list(bound = 12, max_nodes = 1000000L), fixed = c(cor_error = 0))
  expect_equal(z$loglik, unname(log(s['value']*y['value'])), tolerance = 1e-9)
  # The density derivative at rho0 is u*v. Factorization does not zero the
  # free correlation score: its expectation is the product of these moments.
  expect_equal(unname(z$score['cor_re']), unname(s['mean']*y['mean']), tolerance = 1e-9)
  expect_gt(abs(z$score['cor_re']), 1e-5)
})

test_that('physical near-boundary finite values do not bypass accuracy rejection', {
  # The physical Gaussian ridge at rho=.999 is intentionally unresolved at
  # order20. This is a test of rejection, not a claim of uniform GL accuracy.
  data = data.frame(id = 1, d = 1, y = 1)
  prepared = prepare_data(y ~ 1, d ~ 1, data, 'id', 'probit', 'fail')
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0,
    sd_re_selection = 2.2, sd_re_outcome = 9, cor_re = .999, cor_error = 0)
  cache = prepare_split_panels(prepared)
  raw = split_physical_panel(theta, cache$panels[[1]], 20L,
    list(bound = 12, max_nodes = 1000000L, .single_outcome = FALSE), fixed = c(cor_error = 0))
  r = .999*2.2*9/sqrt((1+2.2^2)*(1+9^2))
  expect_gt(abs(raw$loglik-log(.25+asin(r)/(2*pi))), 1e-3)
  expect_error(evaluate_split_binary(theta, prepared,
    panel_control(quadrature = 'gl', nodes = 20), c(cor_error = 0), TRUE, single_outcome = FALSE),
    class = 'panelselect_split_accuracy_error')
  # The default conditional identity eliminates this tensor ridge. Its finite
  # rule is different and is checked against the same exact probability.
  conditional = evaluate_model(theta, prepared,
    panel_control(quadrature = 'gl', nodes = 20), c(cor_error = 0), TRUE)
  expect_equal(conditional$loglik, log(.25 + asin(r)/(2*pi)), tolerance = 1e-9)
  expect_identical(conditional$diagnostics$routes, 'conditional_single_outcome')
})
