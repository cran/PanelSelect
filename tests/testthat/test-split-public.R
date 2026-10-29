test_that('standard split controls are explicit, validated and update-safe', {
  control = panel_control(quadrature = 'gl', nodes = 20)
  expect_identical(control$nodes, 20L)
  expect_identical(control$quadrature_options$reference_nodes, 10L)
  expect_identical(panel_control(quadrature = 'gl')$nodes, 20L)
  expect_identical(panel_control(quadrature = 'gl', nodes = 10)$nodes, 10L)
  expect_identical(validate_panel_control(control), control)
  expect_identical(panel_control()$quadrature, 'gh')
  expect_identical(panel_control()$nodes, 10L)
  expect_error(panel_control(quadrature = 'gl', nodes = c(10, 20)), 'one integer')
  expect_error(panel_control(quadrature = 'gl', nodes = 2), 'one integer')
  expect_error(panel_control(quadrature = 'gl', nodes = 20,
    quadrature_options = list(reference_nodes = 20)), 'different from nodes')
  expect_error(panel_control(quadrature_options = list(bound = 12)), 'require')
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(bound = NA)), 'positive finite')
})

test_that('split fitting uses analytical scores for optimization and native information', {
  skip_if_not_installed('maxLik')
  data = data.frame(id = rep(1:12, each = 2),
    d = rep(c(1,1,1,0,1,1), 4), y = rep(c(1,0,0,NA,1,1), 4))
  theta = c('selection.(Intercept)' = -.2, 'outcome.(Intercept)' = .1,
    sd_re_selection = .8, sd_re_outcome = 1.1, cor_re = .25, cor_error = 0)
  fixed = theta[names(theta) != 'outcome.(Intercept)']
  fit = panel_select(y ~ 1, d ~ 1, data, 'id', 'probit', start = theta, fixed = fixed,
    control = panel_control(quadrature = 'gl', optimizer = 'maxLik', reltol = 1e-10), vcov = 'observed')
  expect_true(fit$optimizer$converged)
  expect_lt(max(abs(fit$optimizer$gradient)), 1e-5)
  expect_true(all(is.finite(vcov(fit))))
  expect_gt(vcov(fit)['outcome.(Intercept)', 'outcome.(Intercept)'], 0)
  accuracy = attr(fit$optimizer$information, 'checked_accuracy')
  expect_identical(accuracy$status, 'PASS')
  expect_match(accuracy$basis, 'Analytical finite moving-cut GL')
  expect_true(all(check_quadrature(fit, nodes = 30)$comparison$within_tolerance))
  expect_identical(fit$integration$gradient_method, 'analytic')
  expect_true(fit$integration$score_checked)
  value = evaluate_model(coef(fit), fit$prepared, fit$control, fit$fixed, score = FALSE)
  expect_false(value$diagnostics$score_checked)
  expect_identical(value$loglik, unname(fit$loglik_by_unit))
})

test_that('standard split fitting has an analytical objective', {
  data = data.frame(id = rep(1:3, each = 2), d = c(1,1,1,0,0,0),
    y = c(1,1,0,NA,NA,NA))
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0,
    sd_re_selection = 1, sd_re_outcome = 1, cor_re = 0, cor_error = 0)
  fit = panel_select(y ~ 1, d ~ 1, data, 'id', 'probit', fixed = theta,
    control = panel_control(quadrature = 'gl', nodes = 20))
  expect_equal(unname(fit$loglik_by_unit), log(c(1/9, 1/12, 1/3)), tolerance = 1e-10)
  expect_identical(fit$integration$method, 'split_legendre')
  expect_identical(fit$integration$gradient_method, 'analytic')
  expect_equal(fit$score, matrix(0, 3, length(theta), dimnames = list(c('1', '2', '3'), names(theta))))
  value = evaluate_model(theta, fit$prepared, fit$control, theta, score = FALSE)
  expect_identical(value$loglik, unname(fit$loglik_by_unit))
  expect_null(value$score)
  expect_output(print(fit), 'Gauss-Legendre \\(gl\\) with moving cuts')
  expect_output(print(summary(fit)), 'Gauss-Legendre \\(gl\\) with moving cuts')
  expect_true(all(check_quadrature(fit, nodes = 30)$comparison$within_tolerance))
  expect_error(check_quadrature(fit, nodes = 10), 'different from nodes|increase')
  expect_error(check_quadrature(fit, observation_nodes = 30), 'no PLN')
  expect_equal(predict(fit, data, conditioning = 'population'),
    setNames(rep(.5, 6), as.character(1:6)), tolerance = 1e-14)
  expect_error(random_effects(fit), 'GL integration')
  expect_error(posterior_rule(theta, fit$prepared, fit$control, theta), 'GL integration')
  updated = update(fit, control = panel_control(quadrature = 'gl', nodes = 30))
  expect_identical(updated$control$quadrature, 'gl')
  expect_equal(updated$loglik, fit$loglik, tolerance = 1e-10)
  free_error = panel_select(y ~ 1, d ~ 1, data, 'id', 'probit',
    start = theta, fixed = theta[names(theta) != 'cor_error'],
    control = panel_control(quadrature = 'gl', nodes = 20, maxit = 0))
  expect_true(is.finite(free_error$loglik))
  expect_true(is.list(free_error$control$integration_mesh))
})

test_that('standalone split probit supports auxiliaries and exact zero-scale models', {
  data = data.frame(id = rep(1:2, each = 2), y = c(1,1,1,0))
  theta = c('outcome.(Intercept)' = 0, sd_re_outcome = 1)
  fit = panel_select(y ~ 1, NULL, data, 'id', 'probit', fixed = theta,
    control = panel_control(quadrature = 'gl', nodes = 20), init = 're')
  expect_equal(unname(fit$loglik_by_unit), log(c(1/3, 1/6)), tolerance = 1e-10)
  independent = panel_select(y ~ 1, NULL, data, 'id', 'probit', fixed = c(theta[1], sd_re_outcome = 0),
    control = panel_control(quadrature = 'gl', nodes = 20))
  expect_equal(unname(independent$loglik_by_unit), rep(log(1/4), 2), tolerance = 1e-14)
  expect_error(check_quadrature(independent), 'no numerical integral')
  poisson = panel_select(y ~ 1, NULL, data, 'id', 'poisson', fixed = c(theta[1], sd_re_outcome = 0),
    control = panel_control(quadrature = 'gl', nodes = 20))
  expect_equal(unname(poisson$loglik_by_unit), c(-2, -2), tolerance = 1e-13)
})

test_that('moving-cut GL budget messages name the rule and the controlling option', {
  data = data.frame(id = rep(1:3, each = 2), d = c(1, 1, 1, 0, 0, 0), y = c(1, 1, 0, NA, NA, NA))
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0, sd_re_selection = 1,
    sd_re_outcome = 1, cor_re = 0, cor_error = 0)
  expect_error(panel_select(y ~ 1, d ~ 1, data, 'id', 'probit', fixed = theta,
    control = panel_control(quadrature = 'gl', nodes = 20, quadrature_options = list(max_nodes = 5))),
    'GL quadrature (moving cuts) node budget exhausted (quadrature_options$max_nodes = 5)', fixed = TRUE)
})
