precision_fit = function(log_sd = -.7, correlation_eta = .3) {
  info = data.frame(name = c('outcome.(Intercept)', 'sd_re_outcome', 'cor_re'),
    type = c('coefficient', 'sd', 'correlation'), fixed = FALSE)
  eta = stats::setNames(c(.2, log_sd, correlation_eta), info$name)
  theta = parameter_from_optimizer(eta, eta, info)
  V = diag(c(.04, 1, .09)); dimnames(V) = list(info$name, info$name)
  J = parameter_jacobian(theta, info)
  natural = V*outer(J, J)
  structure(list(coefficients = theta, parameter_info = info,
    covariance = list(observed = natural, opg = 2*natural, sandwich = 3*natural),
    optimizer = list(par = eta, covariance = list(observed = V, opg = 2*V, sandwich = 3*V),
      information_computed = TRUE, hessian = -diag(3), gradient = eta*0,
      counts = c('function' = 1L), code = 0L, converged = TRUE),
    control = list(vcov = 'observed'), fixed = numeric(),
    score = matrix(0, 2, 3), loglik = -1, prepared = list(n_rows = 4L, n_panels = 2L)),
    class = 'panelselect')
}

test_that('retained coordinates preserve tiny SD SEs and intervals', {
  fit = precision_fit(-400, 18)
  answer = coefficient_table(fit, 'observed')
  expect_identical(vcov(fit), fit$covariance$observed)
  expect_equal(vcov(fit)['sd_re_outcome', 'sd_re_outcome'], 0)
  expect_gt(answer['sd_re_outcome', 'se'], 0)
  expect_equal(answer['sd_re_outcome', 'se']/exp(-400), 1)
  expect_equal(unname(answer['sd_re_outcome', c('lci','uci')]/
    exp(-400+c(-1,1)*qnorm(.975))), c(1,1))
  expect_true(is.na(answer['sd_re_outcome','z']))
  expect_equal(unname(answer['cor_re','z']), 18/.3, tolerance = 1e-14)
  expect_gt(abs(atanh(fit$coefficients['cor_re'])-18), .001)
  expect_equal(confint(fit, 'sd_re_outcome')[1, ]/exp(-400+c(-1,1)*qnorm(.975)),
    stats::setNames(c(1,1), colnames(confint(fit, 'sd_re_outcome'))))
})

test_that('named retained scales agree ordinarily and reject stale metadata', {
  fit = precision_fit(); retained = coefficient_table(fit, 'observed')
  old = fit; old$optimizer$covariance = NULL; old$optimizer$par = NULL
  expect_equal(coefficient_table(old, 'observed'), retained, tolerance = 1e-14)
  shuffled = fit; index = c(3,1,2)
  shuffled$optimizer$par = shuffled$optimizer$par[index]
  shuffled$optimizer$covariance$observed = shuffled$optimizer$covariance$observed[index,rev(index)]
  expect_equal(coefficient_table(shuffled, 'observed'), retained)
  stale = fit; stale$optimizer$par[1] = stale$optimizer$par[1]+1
  expect_identical(coefficient_table(stale, 'observed'), coefficient_table(old, 'observed'))
  stale = fit; dimnames(stale$optimizer$covariance$observed) = NULL
  expect_identical(coefficient_table(stale, 'observed'), coefficient_table(old, 'observed'))
  stale = fit; names(stale$optimizer$par)[2] = names(stale$optimizer$par)[1]
  expect_identical(coefficient_table(stale, 'observed'), coefficient_table(old, 'observed'))
  stale = fit; stale$optimizer$par = stale$optimizer$par+0i
  expect_identical(coefficient_table(stale, 'observed'), coefficient_table(old, 'observed'))
  stale = fit; stale$optimizer$covariance$observed = stale$optimizer$covariance$observed+0i
  expect_identical(coefficient_table(stale, 'observed'), coefficient_table(old, 'observed'))
  permuted = old; permuted$covariance$observed = permuted$covariance$observed[index,rev(index)]
  expect_identical(coefficient_table(permuted, 'observed'), coefficient_table(old, 'observed'))
  expect_identical(vcov(permuted), permuted$covariance$observed)
  incompatible = old; rownames(incompatible$covariance$observed)[1] = 'wrong'
  expect_error(coefficient_table(incompatible, 'observed'), 'names matching coefficients')
  incompatible = old; dimnames(incompatible$covariance$observed) = NULL
  expect_error(coefficient_table(incompatible, 'observed'), 'names matching coefficients')
  incompatible = old; rownames(incompatible$covariance$observed)[2] = rownames(incompatible$covariance$observed)[1]
  expect_error(coefficient_table(incompatible, 'observed'), 'names matching coefficients')
  incompatible = old; incompatible$covariance$observed = incompatible$covariance$observed+0i
  expect_error(coefficient_table(incompatible, 'observed'), 'must be real')
})

test_that('fixed and unavailable inference remains unavailable without substitution', {
  fit = precision_fit()
  fit$covariance$observed[] = NA_real_; fit$optimizer$covariance$observed[] = NA_real_
  local_mocked_bindings(information_matrix = function(...) stop('unexpected Hessian'), .package = 'PanelSelect')
  expect_true(all(is.na(coefficient_table(fit, 'observed')[, 'se'])))
  expect_true(all(is.finite(coefficient_table(fit, 'opg')[, 'se'])))
  fit = precision_fit(); fit$parameter_info$fixed[] = TRUE
  fit$coefficients['sd_re_outcome'] = 0
  fit$optimizer$par = numeric(); fit$optimizer$covariance$observed = matrix(numeric(),0,0)
  fit$covariance$observed[] = 0
  answer = coefficient_table(fit, 'observed')
  expect_true(all(is.na(answer[, c('se','z','p','lci','uci')])))
  expect_identical(unname(answer[, 'estimate']), unname(fit$coefficients))
})

test_that('on-demand observed and sandwich tables compute one Hessian each', {
  fit = precision_fit(-400); saved = fit
  fit$control$vcov = 'opg'; fit$optimizer$information_computed = FALSE
  fit$optimizer$covariance[c('observed','sandwich')] = list(NULL,NULL)
  fit$covariance$observed[] = NA_real_; fit$covariance$sandwich[] = NA_real_
  calls = 0L
  local_mocked_bindings(fit_evaluator = function(...) list(),
    information_matrix = function(...) { calls <<- calls+1L; diag(3) },
    covariance_matrices = function(...) list(natural = saved$covariance,
      optimizer = saved$optimizer$covariance), .package = 'PanelSelect')
  expect_gt(coefficient_table(fit, 'observed')['sd_re_outcome','se'], 0)
  expect_identical(calls, 1L)
  expect_gt(coefficient_table(fit, 'sandwich')['sd_re_outcome','se'], 0)
  expect_identical(calls, 2L)
  coefficient_table(fit, 'opg'); expect_identical(calls, 2L)
  expect_false(fit$optimizer$information_computed)
})
