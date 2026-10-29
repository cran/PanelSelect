# Public method options: fitted(), confint(level =),
# AIC(k =), random uniform starts and other optimizer methods.

methods_fixture = function() {
  local_test_seed(1601)
  data = data.frame(id = rep(1:20, each = 3), x = stats::rnorm(60))
  data$y = stats::rpois(60, exp(.2 + .4*data$x + .5*stats::rnorm(20)[data$id]))
  data
}

test_that('fitted values, confidence levels and AIC penalties follow their definitions', {
  data = methods_fixture()
  fit = panel_select(y ~ x, NULL, data, 'id', 'poisson', control = panel_control(nodes = 6L))
  expect_identical(fitted(fit), predict(fit, type = 'response'))
  table = summary(fit)$coefficients
  critical = stats::qnorm(.95)
  estimate = table[, 'estimate']; se = table[, 'se']
  expected = cbind(estimate - critical*se, estimate + critical*se)
  # Scales use log-scale intervals: their optimizer-scale SE is se/estimate.
  expected['sd_re_outcome', ] = exp(log(estimate[['sd_re_outcome']]) +
    c(-1, 1)*critical*se[['sd_re_outcome']]/estimate[['sd_re_outcome']])
  interval = confint(fit, level = .9)
  expect_identical(colnames(interval), c('5 %', '95 %'))
  expect_equal(unname(interval), unname(expected), tolerance = 1e-12)
  df = attr(logLik(fit), 'df')
  expect_equal(AIC(fit, k = 3), -2*as.numeric(logLik(fit)) + 3*df, tolerance = 1e-14)
})

test_that('random uniform starts and optim methods reach the BFGS optimum', {
  data = methods_fixture()
  # Hold the finite AGH objective fixed across starts and optimizer methods.
  # A tightly converged BFGS fit is the reference optimum.
  reference = panel_select(y ~ x, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'agh', nodes = 6L, reltol = 1e-12))
  local_test_seed(1602)
  uniform = panel_select(y ~ x, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'agh', nodes = 6L), init = 'unif')
  expect_true(uniform$optimizer$converged)
  expect_equal(coef(uniform), coef(reference), tolerance = 1e-4)
  for (method in c('CG', 'L-BFGS-B')) {
    fit = panel_select(y ~ x, NULL, data, 'id', 'poisson',
      start = coef(reference) + c(.05, -.05, .05),
      control = panel_control(quadrature = 'agh', nodes = 6L, method = method, maxit = 2000L))
    expect_true(fit$optimizer$converged, info = method)
    expect_equal(coef(fit), coef(reference), tolerance = 1e-4, info = method)
  }
})

test_that('native maxLik Newton reaches the BFGS optimum', {
  skip_if_not_installed('maxLik')
  data = methods_fixture()
  reference = panel_select(y ~ x, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'agh', nodes = 6L, reltol = 1e-12))
  newton = panel_select(y ~ x, NULL, data, 'id', 'poisson', start = coef(reference) + c(.05, -.05, .05),
    control = panel_control(quadrature = 'agh', nodes = 6L, optimizer = 'maxLik', method = 'NR'))
  expect_true(newton$optimizer$converged)
  expect_equal(coef(newton), coef(reference), tolerance = 1e-4)
})
