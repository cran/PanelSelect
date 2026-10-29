for (backend in c('optim', 'maxLik')) {
  test_that(paste('default inference is independent-unit OPG without information:', backend), {
    if (backend == 'maxLik') skip_if_not_installed('maxLik')
    expect_identical(eval(formals(panel_select)$vcov)[[1L]], 'opg')
    data = data.frame(id = rep(seq_len(6), each = 2),
      y = c(1,2,0,3,1,0,2,1,4,0,1,3),
      x = c(-.8,.2,-.4,.7,.1,-.2,.5,-.3,.9,-.6,.3,.4))
    calls = 0L; original = information_matrix
    local_mocked_bindings(information_matrix = function(...) {
      calls <<- calls + 1L; original(...)
    }, .package = 'PanelSelect')
    fit = panel_select(y ~ x, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0),
      control = panel_control(optimizer = backend))
    expect_identical(fit$control$vcov, 'opg')
    expect_false(fit$optimizer$information_computed)
    expect_null(fit$optimizer$information_error)
    expect_null(fit$optimizer$hessian)
    expect_true(all(is.na(fit$optimizer$information)))
    free = !fit$parameter_info$fixed
    expected = solve(crossprod(fit$score[, free, drop = FALSE]))
    expect_equal(vcov(fit)[free, free, drop = FALSE], expected, tolerance = 1e-12)
    before = calls
    expect_identical(check_stationarity(fit)$status, 'unavailable')
    expect_equal(summary(fit)$coefficients[, 'se'][free], sqrt(diag(expected)), tolerance = 1e-12)
    expect_true(all(is.finite(confint(fit, parm = names(coef(fit))[free]))))
    expect_true(all(is.finite(predict(fit, data, se.fit = TRUE)$se.fit)))
    expect_identical(calls, before)
    expect_true(all(is.finite(vcov(fit, type = 'observed'))))
    expect_identical(calls, before + 1L)
    expect_true(all(is.finite(vcov(fit, type = 'sandwich'))))
    expect_identical(calls, before + 2L)
    observed = panel_select(y ~ x, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0),
      control = panel_control(optimizer = backend), vcov = 'observed')
    expect_identical(calls, before + 3L)
    expect_true(observed$optimizer$information_computed)
    expect_identical(coef(fit), coef(observed))
    expect_identical(fit$loglik, observed$loglik)
  })
}

test_that('OPG covariance transforms panel scores without observed information', {
  theta = c(b = .4, sd = 2, rho = .3, fixed = 0)
  info = data.frame(name = names(theta), type = c('coefficient','sd','correlation','coefficient'),
                    fixed = c(FALSE,FALSE,FALSE,TRUE))
  score = matrix(c(1,2,-1,0, -.4,1,2,0, 2,-1,.5,0, .2,.3,1.1,0), 4, byrow = TRUE,
                 dimnames = list(NULL, names(theta)))
  jacobian = c(1,2,1-.3^2)
  eta_score = sweep(score[, 1:3], 2, jacobian, '*')
  expected = solve(crossprod(eta_score))
  covariance = covariance_matrices(theta, info, score,
    matrix(NA_real_, 3, 3, dimnames = list(names(theta)[1:3], names(theta)[1:3])))
  expect_equal(covariance$optimizer$opg, expected, tolerance = 1e-12)
  expect_equal(covariance$natural$opg[1:3,1:3], expected * outer(jacobian,jacobian), tolerance = 1e-12)
  expect_true(all(covariance$natural$opg[4,] == 0))
  expect_true(all(is.na(covariance$optimizer$observed)))
  expect_true(all(is.na(covariance$optimizer$sandwich)))
})

test_that('all maxLik OPG paths suppress the final Hessian', {
  skip_if_not_installed('maxLik')
  calls = list()
  local_mocked_bindings(maxLik = function(logLik, grad, start, method, control, ...) {
    calls[[method]] <<- list(...)
    list(estimate = start, code = 0L, iterations = 1L, message = 'test')
  }, .package = 'maxLik')
  evaluator = list(evaluate = function(x) list(loglik = -x^2,
    optimizer_score = matrix(-2*x, 1L)))
  for (method in c('BFGS','NR','BHHH','NM','CG')) {
    optimize_model(.2, evaluator, fit_control(panel_control(optimizer = 'maxLik', method = method)))
    expect_identical(calls[[method]], list(finalHessian = FALSE))
  }
})
