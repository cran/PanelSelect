test_that('higher split references check rather than replace the canonical rule', {
  data = data.frame(id = rep(1:3, each = 2), d = c(1,1,1,0,0,0), y = c(1,0,0,NA,NA,NA))
  theta = c('selection.(Intercept)' = -.2, 'outcome.(Intercept)' = .1,
    sd_re_selection = .8, sd_re_outcome = 1.1, cor_re = .25, cor_error = 0)
  prepared = prepare_data(y ~ 1, d ~ 1, data, 'id', 'probit', 'fail')
  lower = evaluate_model(theta, prepared, panel_control(quadrature = 'gl', nodes = 12,
    quadrature_options = list(reference_nodes = 10)), c(cor_error = 0))
  higher = evaluate_model(theta, prepared, panel_control(quadrature = 'gl', nodes = 12,
    quadrature_options = list(reference_nodes = 15)), c(cor_error = 0))
  expect_identical(higher$loglik, lower$loglik)
  expect_identical(higher$score, lower$score)
  expect_identical(higher$diagnostics$canonical_order, 12L)
  expect_identical(higher$diagnostics$reference_order, 15L)
  expect_error(panel_control(quadrature = 'gl', nodes = 12,
    quadrature_options = list(reference_nodes = 257)), '2 through 256')
})

test_that('split Linear probabilities and analytical scores match independent integrals', {
  data = data.frame(id = rep(1:3, each = 3), x = rep(c(-.4, .1, .6), 3),
    d = c(1,1,1,1,0,1,0,0,0), y = c(.2,-.7,1.1,.8,NA,-.2,NA,NA,NA))
  theta = c('selection.(Intercept)' = -.3, 'selection.x' = .5,
    'outcome.(Intercept)' = .2, 'outcome.x' = -.4,
    sd_re_selection = 1.3, sd_re_outcome = 1.7, sd_error = .9, cor_re = .4, cor_error = 0)
  prepared = prepare_data(y ~ x, d ~ x, data, 'id', 'gaussian', 'fail')
  control = panel_control(quadrature = 'gl', nodes = 12, quadrature_options = list(reference_nodes = 15))
  # Original two-dimensional density; no outcome-conditioning implementation.
  original = function(rows, par) {
    selected = rows[data$d[rows] == 1]
    a = par[1] + par[2]*data$x[rows]
    b = par[3] + par[4]*data$x[selected]
    conditional = function(v) vapply(v, function(vv) {
      selection = integrate(function(u) {
        prob = vapply(u, function(uu) prod(pnorm((2*data$d[rows]-1)*(a+par['sd_re_selection']*uu))), numeric(1))
        prob*dnorm(u, par['cor_re']*vv, sqrt(1-par['cor_re']^2))
      }, -10, 10, rel.tol = 1e-10, abs.tol = 1e-12)$value
      selection*prod(dnorm(data$y[selected], b+par['sd_re_outcome']*vv, par['sd_error']))*dnorm(vv)
    }, numeric(1))
    log(integrate(conditional, -10, 10, rel.tol = 1e-9, abs.tol = 1e-12)$value)
  }
  for (rho in c(-.6, 0, .4)) {
    theta['cor_re'] = rho
    value = evaluate_model(theta, prepared, control, c(cor_error = 0))
    reference = vapply(split(seq_len(nrow(data)), data$id), original, numeric(1), par = theta)
    expect_equal(value$loglik, reference, tolerance = 2e-8, ignore_attr = TRUE)
    for (j in seq_along(theta)[names(theta) != 'cor_error']) {
      step = 1e-5*max(1, abs(theta[j])); plus = minus = theta
      plus[j] = plus[j]+step; minus[j] = minus[j]-step
      numerical = (evaluate_model(plus, prepared, control, c(cor_error = 0), FALSE)$loglik -
        evaluate_model(minus, prepared, control, c(cor_error = 0), FALSE)$loglik)/(2*step)
      expect_equal(value$score[, j], numerical, tolerance = 2e-6, ignore_attr = TRUE)
    }
    if (rho == 0) expect_gt(abs(sum(value$score[, 'cor_re'])), .01)
    # Panel history deduplication preserves the original order and marginals.
    subset = subset_panels(prepared, c(3L, 1L, 2L))
    reordered = evaluate_model(theta, subset, control, c(cor_error = 0))
    expect_equal(reordered$loglik, value$loglik[c(3,1,2)], tolerance = 1e-12)
    expect_equal(reordered$score, value$score[c(3,1,2), ], tolerance = 1e-12)
  }
})

test_that('split Linear one-wave likelihood uses the exact conditional-normal identity', {
  data = data.frame(id = 1:3, d = c(1,1,0), y = c(.8, -.7, NA))
  theta = c('selection.(Intercept)' = -.4, 'outcome.(Intercept)' = .2,
    sd_re_selection = 1.2, sd_re_outcome = 2, sd_error = .8, cor_re = -.5, cor_error = 0)
  prepared = prepare_data(y ~ 1, d ~ 1, data, 'id', 'gaussian', 'fail')
  control = panel_control(quadrature = 'gl', nodes = 12, quadrature_options = list(reference_nodes = 15))
  value = evaluate_model(theta, prepared, control, c(cor_error = 0))
  variance = theta['sd_re_outcome']^2+theta['sd_error']^2
  mean_u = theta['cor_re']*theta['sd_re_outcome']*(data$y[1:2]-theta[2])/variance
  variance_u = 1-theta['cor_re']^2*theta['sd_re_outcome']^2/variance
  reference = c(dnorm(data$y[1:2], theta[2], sqrt(variance), log = TRUE) +
    pnorm((theta[1]+theta['sd_re_selection']*mean_u)/sqrt(1+theta['sd_re_selection']^2*variance_u), log.p = TRUE),
    pnorm(-theta[1]/sqrt(1+theta['sd_re_selection']^2), log.p = TRUE))
  expect_equal(value$loglik, unname(reference), tolerance = 1e-10)
  free_tau = panel_select(y ~ 1, d ~ 1, data, 'id', 'gaussian', start = theta,
    fixed = theta[names(theta) != 'cor_error'],
    control = panel_control(quadrature = 'gl', nodes = 12, maxit = 0,
      quadrature_options = list(reference_nodes = 15)))
  expect_true(all(is.finite(free_tau$loglik_by_unit)))
  zero = theta; zero['sd_re_selection'] = 0
  independent = evaluate_model(zero, prepared, control, c(cor_error = 0, sd_re_selection = 0))
  expect_true(all(is.finite(independent$loglik)))
  zero = theta; zero['sd_re_outcome'] = 0
  independent = evaluate_model(zero, prepared, control, c(cor_error = 0, sd_re_outcome = 0))
  expect_true(all(is.finite(independent$loglik)))
  tiny = theta; tiny['sd_error'] = 1e-200
  expect_error(evaluate_model(tiny, prepared, control, c(cor_error = 0)),
    class = 'panelselect_integration_error')
})

test_that('native split Linear fitting and population uncertainty use OPG', {
  skip_if_not_installed('maxLik')
  data = data.frame(id = rep(1:12, each = 2), d = rep(c(1,1,1,0), 6),
    y = rep(c(.8,-.4,.2,NA), 6))
  theta = c('selection.(Intercept)' = -.2, 'outcome.(Intercept)' = .1,
    sd_re_selection = .8, sd_re_outcome = 1.1, sd_error = .7, cor_re = .25, cor_error = 0)
  fixed = theta[names(theta) != 'outcome.(Intercept)']
  fit = panel_select(y ~ 1, d ~ 1, data, 'id', 'gaussian', start = theta, fixed = fixed,
    control = panel_control(quadrature = 'gl', nodes = 12,
      optimizer = 'maxLik', reltol = 1e-10, quadrature_options = list(reference_nodes = 15)), init = 're')
  expect_true(fit$optimizer$converged)
  expect_false(fit$optimizer$information_computed)
  expect_equal(vcov(fit), fit$covariance$opg)
  expect_true(all(is.finite(predict(fit, data, conditioning = 'population', se.fit = TRUE)$se.fit)))
  expect_true(all(check_quadrature(fit, nodes = 15)$comparison$within_tolerance))
})
