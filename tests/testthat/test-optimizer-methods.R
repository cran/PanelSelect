test_that('optim methods are validated and stored by their full names', {
  expect_identical(panel_control(method = 'Nelder')$method, 'Nelder-Mead')
  expect_identical(panel_control(method = 'L-BFGS')$method, 'L-BFGS-B')
  expect_identical(panel_control(method = 'BFGS')$method, 'BFGS')
  expect_error(panel_control(method = 'foo'), 'method must be one of')
  expect_error(panel_control(method = 'Brent'), 'method must be one of')
  # Exact-name option rules now apply to an abbreviated L-BFGS-B request.
  options = optimizer_options(panel_control(method = 'L-BFGS'))
  expect_true('factr' %in% names(options))
  expect_false('reltol' %in% names(options))
})

test_that('maxLik methods are validated without changing accepted spellings', {
  for (method in c('BFGS', 'bfgs', 'NR', 'Newton-Raphson', 'BHHH', 'BFGS-R', 'BFGSR',
                   'CG', 'Conjugate-Gradient', 'NM', 'Nelder-Mead', 'SANN'))
    expect_identical(panel_control(optimizer = 'maxLik', method = method)$method, method)
  expect_error(panel_control(optimizer = 'maxLik', method = 'foo'), 'method must be one of')
})

test_that('simulated annealing is not reported as converged', {
  # optim and maxLik return code 0 for SANN whatever the stopping state.
  for (backend in c('optim', 'maxLik')) {
    control = panel_control(optimizer = backend, method = 'SANN')
    options = optimizer_options(control)
    for (code in c(0, 1)) expect_identical(optimizer_status(code, control, options, 1L),
      'no_convergence_test')
  }
})

test_that('maxLik iteration counts carry plain names', {
  skip_if_not_installed('maxLik')
  data = data.frame(id = rep(1:6, each = 3L), x = sin(1:18),
    d = c(1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1))
  fit = panel_select(d ~ x, NULL, data, 'id', 'probit',
    control = panel_control(quadrature = 'gh', nodes = 5L, optimizer = 'maxLik', maxit = 0L))
  expect_identical(names(fit$optimizer$counts), c('function', 'gradient'))
})
