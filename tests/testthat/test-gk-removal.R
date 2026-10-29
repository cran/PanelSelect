# Unsupported GK settings stop with a pointer to GL. Historical fits and a
# check in fixtures/saved-gk.R still display, stop where they would integrate
# with GK, and refit with GL.

gk = dget(test_path('fixtures', 'saved-gk.R'))

test_that('the removed GK rule, its option and its check stop with a pointer to GL', {
  expect_identical(eval(formals(panel_control)$quadrature), c('gh', 'agh', 'gl'))
  expect_identical(eval(formals(panel_control)$observation_quadrature), c('gh', 'agh', 'gl'))
  expect_identical(eval(formals(check_quadrature)$compare), c('higher_order', 'gl'))
  expect_error(panel_control(quadrature = 'gk'), 'quadrature = "gk" was removed; use "gl"',
    fixed = TRUE)
  expect_error(panel_control(observation_quadrature = 'gk'),
    'observation_quadrature = "gk" was removed; use "gl"', fixed = TRUE)
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(max_bound = 32)),
    'quadrature_options entry max_bound belonged to the removed GK rule; GL keeps its bound',
    fixed = TRUE)
  fit = panel_select(count ~ x, NULL, gk$data$count, 'id', 'poisson',
    control = panel_control(quadrature = 'gh'))
  expect_error(check_quadrature(fit, compare = 'gk'),
    'compare = "gk" was removed; use compare = "gl"', fixed = TRUE)
})

test_that('fits saved with GK display, stop where they would integrate, and refit with GL', {
  count = gk$fits$count_gk
  hint = 'refit with update(fit, control = panel_control(quadrature = "gl"))'
  expect_identical(names(coef(count)), c('outcome.(Intercept)', 'outcome.x', 'sd_re_outcome'))
  expect_true(all(is.finite(vcov(count))))
  expect_no_error(utils::capture.output(print(count)))
  printed = utils::capture.output(print(summary(count)))
  expect_true(any(grepl('Gauss-Kronrod 21/10 (gk, removed)', printed, fixed = TRUE)))
  expect_error(update(count), hint, fixed = TRUE)
  expect_error(check_quadrature(count), hint, fixed = TRUE)
  expect_error(add_information(count), hint, fixed = TRUE)
  refit = update(count, control = panel_control(quadrature = 'gl'), init = 'pooled')
  expect_identical(refit$control$quadrature, 'gl')
  expect_identical(refit$optimizer$status, 'converged')
  # Both rules meet the staged total log-likelihood tolerance (1e-5) at their
  # optima, which therefore agree to within twice that.
  expect_lt(abs(refit$loglik - count$loglik), 2e-5)

  pln = gk$fits$pln_gk
  both = 'refit with update(fit, control = panel_control(quadrature = "gl", observation_quadrature = "gl"))'
  expect_no_error(utils::capture.output(print(summary(pln))))
  expect_error(update(pln), both, fixed = TRUE)
  expect_error(check_quadrature(pln), both, fixed = TRUE)
  # A check that compared with GK still prints.
  expect_match(utils::capture.output(print(gk$check))[1L], '^Quadrature check [(]gk[)] of: ')
})

# compare = "gl" replaces compare = "gk"; its agreement with the frozen refined
# references is tested in test-interface-capabilities.R.
test_that('compare = "gl" needs a fit without GL', {
  fit = panel_select(count ~ x, NULL, gk$data$count, 'id', 'poisson',
    control = panel_control(quadrature = 'gl'))
  expect_error(check_quadrature(fit, compare = 'gl'),
    'this fit already uses GL; use compare = "higher_order"', fixed = TRUE)
})

test_that('GL keeps its bound: a far-tail panel stops with the bound hint and a larger bound resolves it', {
  # Counts near 200 with a zero intercept and person-effect SD 0.4 put the
  # posterior mode of the standard-normal person coordinate at 13.1, beyond
  # the default bound 12, which GK widened on demand.
  data = data.frame(id = 1L, y = c(200, 210, 190))
  theta = c('outcome.(Intercept)' = 0, sd_re_outcome = .4)
  fit = function(bound) panel_select(y ~ 1, NULL, data, 'id', 'poisson', fixed = theta,
    control = panel_control(quadrature = 'gl', quadrature_options = list(bound = bound)))
  expect_error(fit(12), paste('Outer omitted-tail envelope unresolved at maximum bound',
    '(maximum bound 12; raise quadrature_options$bound)'), fixed = TRUE)
  # Independent reference: the panel likelihood integrated over its posterior
  # region, scaled near its peak to avoid underflow.
  log_f = function(z) vapply(z, function(v) sum(stats::dpois(data$y, exp(.4*v), log = TRUE)) +
    stats::dnorm(v, log = TRUE), numeric(1))
  peak = log_f(13.1)
  reference = peak + log(stats::integrate(function(z) exp(log_f(z) - peak), 9, 17, rel.tol = 1e-12)$value)
  # The default per-panel log-likelihood tolerance is 1e-7.
  expect_lt(abs(fit(16)$loglik - reference), 1e-7)
})

test_that('the staged time floor is the GL default', {
  control = fit_control(panel_control(quadrature = 'gh', observation_quadrature = 'gl',
    observation_options = list(seconds_total = 10)))
  expect_identical(gl_options(20L)$seconds_total, 300)
  expect_identical(staged_total_seconds(control, staged_options(control, TRUE)), 300)
})
