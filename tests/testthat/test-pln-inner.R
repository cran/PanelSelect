test_that('inner normal derivatives retain precision at the tail switch', {
  x = c(-4.000001, -4, -3.999999, -8, -20, -100, 0, 2)
  # Independent mpmath 1.3.0 derivatives of log(erfc(-x/sqrt(2))/2), 65 digits.
  fourth = c(.009506570178825397, .009506576431595869, .009506582684371073,
    .001105292095455163, .00003570355158845659, 5.988020962773754e-8,
    .1147706820542189, -.1878546856116097)
  fifth = c(.006252768106942366, .006252772838129979, .006252777569321483,
    .0004820640164648541, 6.968588151584452e-6, 2.392816762787677e-9,
    -.004437688462617821, -.03109190219518375)
  value = pln_inner_normal(x, 5L)
  expect_equal(value[[5L]], fourth, tolerance = 2e-9)
  expect_equal(value[[6L]], fifth, tolerance = 2e-9)
  for (order in 1:5) {
    medium = pln_inner_normal(c(-4.000001, -8, -30, -100), 5L, 64L)[[order+1L]]
    deep = pln_inner_normal(c(-4.000001, -8, -30, -100), 5L, 96L)[[order+1L]]
    expect_equal(medium, deep, tolerance = 2e-13)
  }
})

test_that('adaptive inner mixed derivatives match independent high precision', {
  x = local_pln_adaptive(-.3, 4, .8, 3L, lgamma(5), normal_rule(5L, 1L),
                          a = .2, tau = -.4, mode_tol = 1e-12)
  # Independent finite-rule mixed derivatives; the recipe is in fixtures/README.md.
  expect_equal(x$hessian[1, 1, 2], -.2155729636950900, tolerance = 2e-11)
  expect_equal(x$third[1, 1, 1, 3], -.09120888591231401, tolerance = 2e-10)
  expect_equal(x$third[1, 1, 2, 4], .6758976368282991, tolerance = 2e-10)
  expect_equal(x$third[1, 2, 2, 3], .7481695605980433, tolerance = 2e-10)
  expect_equal(x$third[1, 1, 2, 1], .04542390200999184, tolerance = 2e-10)
  expect_equal(x$third[1, 2, 1, 4], .6758976368282991, tolerance = 2e-10)
  expect_gte(x$inner_diagnostics$min_curvature, 1)
  expect_lt(x$inner_diagnostics$max_scaled_residual, 1e-12)
})

test_that('adaptive inner scores differentiate the moving finite objective', {
  evaluate = function(theta, y, order) local_pln_adaptive(theta[2], y, theta[3], order,
    lgamma(y+1), normal_rule(7L, 1L), a = theta[1], tau = theta[4], mode_tol = 1e-12)
  cases = list(c(.3, -.2, .7, .45), c(-.8, .5, 1.1, -.55), c(.2, -.5, .4, 0))
  for (i in seq_along(cases)) {
    theta = cases[[i]]; y = c(3, 9, 0)[i]
    value = evaluate(theta, y, 3L)
    for (step in c(1e-4, 1e-5)) for (j in seq_along(theta)) {
      plus = minus = theta; plus[j] = plus[j]+step; minus[j] = minus[j]-step
      score = (evaluate(plus, y, 0L)$value-evaluate(minus, y, 0L)$value)/(2*step)
      expect_equal(value$gradient[1, j], score, tolerance = 3e-7)
      third = (evaluate(plus, y, 2L)$hessian[1, , 1:2]-
        evaluate(minus, y, 2L)$hessian[1, , 1:2])/(2*step)
      expect_equal(value$third[1, , , j], third, tolerance = 5e-6)
    }
  }
})

test_that('inner batches and requested orders preserve values and axes', {
  a = c(-.5, .2, .7); b = c(-.3, 1, -1); y = c(0, 4, 9)
  rule = normal_rule(9L, 1L)
  batch = local_pln_adaptive(b, y, .8, 3L, lgamma(y+1), rule, a, -.4)
  for (order in 0:3) {
    value = local_pln_adaptive(b, y, .8, order, lgamma(y+1), rule, a, -.4)
    expect_equal(value$value, batch$value, tolerance = 2e-10)
    if (order) expect_equal(value$gradient, batch$gradient, tolerance = 2e-10)
    if (order >= 2L) expect_equal(value$hessian, batch$hessian, tolerance = 2e-9)
    for (i in seq_along(y)) {
      one = local_pln_adaptive(b[i], y[i], .8, order, lgamma(y[i]+1), rule, a[i], -.4)
      expect_equal(value$value[i], one$value, tolerance = 2e-10)
      if (order) expect_equal(value$gradient[i, ], one$gradient[1, ], tolerance = 2e-10)
    }
  }
  standalone = local_pln_adaptive(b, y, .8, 3L, lgamma(y+1), rule)
  selected = local_pln_adaptive(b, y, .8, 3L, lgamma(y+1), rule, a, 0)
  expect_equal(selected$value, standalone$value+pnorm(a, log.p = TRUE), tolerance = 2e-10)
  expect_equal(selected$gradient[, 2:3], standalone$gradient, tolerance = 2e-10)
  expect_true(any(abs(selected$gradient[, 4]) > 1e-5))
  expect_identical(dim(standalone$third), c(3L, 1L, 1L, 2L))
})

test_that('adaptive inner failures remain explicit', {
  rule = normal_rule(7L, 1L)
  expect_error(local_pln_adaptive(-.3, 4, .8, 1L, lgamma(5), rule, .2, -.4,
    mode_maxit = 1L), class = 'panelselect_integration_error')
  expect_error(local_pln_adaptive(-.3, 4, .8, 1L, lgamma(5), rule, .2, 1),
    class = 'panelselect_integration_error')
})
