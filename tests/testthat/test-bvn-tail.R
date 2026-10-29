# High-precision references for log P(X < a, Y < b); see fixtures/README.md.
# kappa is the
# condition number of log P in the double-precision inputs.
bvn_tail_fixture = dget(testthat::test_path('fixtures', 'bvn-tail.R'))$points

# Both batched tail integrals use the Gauss-Kronrod 21/10 rule.
test_that('K21 and embedded G10 have positive rules and their known moments', {
  rule = gk21_rule()
  expect_true(all(rule$weights > 0))
  expect_true(all(rule$reference_weights >= 0))
  expect_equal(sum(rule$weights), 2, tolerance = 2e-15)
  expect_equal(sum(rule$reference_weights), 2, tolerance = 2e-15)
  expect_equal(sum(rule$reference_weights > 0), 10L)
  for (degree in 0:31) {
    truth = if (degree %% 2L) 0 else 2/(degree+1L)
    expect_equal(sum(rule$weights*rule$nodes^degree), truth, tolerance = 4e-15)
    if (degree <= 19L)
      expect_equal(sum(rule$reference_weights*rule$nodes^degree), truth, tolerance = 4e-15)
  }
  expect_equal(rule$nodes, -rev(rule$nodes), tolerance = 0)
})

test_that('the vectorized tail agrees with high-precision references', {
  p = bvn_tail_fixture
  tail = p$region != 'threshold_above'
  value = log_bvn_value(p$a, p$b, p$r)
  # Error allowance: the integration target, or input conditioning in double.
  bound = pmax(2e-11, 4*.Machine$double.eps*p$kappa)
  expect_true(all(abs(value[tail] - p$log_p[tail]) <= bound[tail]))
  # Elements above the 1e-9 threshold keep the pbivnorm path bitwise.
  above = which(!tail)
  expect_identical(value[above], log(pbivnorm::pbivnorm(p$a[above], p$b[above], p$r[above])))
})

test_that('the vectorized tail agrees with the scalar tail routine', {
  grid = expand.grid(a = c(-40, -25, -12, -8, -6), b = c(-30, -10, -6, 0, 4),
    r = c(-.9999, -.95, -.5, -.1, .1, .5, .95, .9999))
  grid = grid[pbivnorm::pbivnorm(grid$a, grid$b, grid$r) < 1e-9, ]
  batch = log_bvn_tail_batch(grid$a, grid$b, grid$r)
  expect_lte(mean(is.na(batch)), .05)
  scalar = mapply(log_bvn_tail, grid$a, grid$b, grid$r)
  # Each path meets a 2e-11 relative integration tolerance.
  expect_lt(max(abs(batch - scalar), na.rm = TRUE), 4e-11)
})

test_that('elements the batch declines use the unchanged scalar routine', {
  # Extreme arguments are left to the scalar routine and its failure reports.
  expect_true(is.na(log_bvn_tail_batch(-2e6, -1, .5)))
  expect_identical(log_bvn_value(-2e6, -1, .5), log_bvn_tail(-2e6, -1, .5))
  values = log_bvn_value(c(-8, -2e6, -9), c(-7, -1, 2), c(-.3, .5, .8))
  expect_identical(values[2L], log_bvn_tail(-2e6, -1, .5))
  expect_equal(values[c(1L, 3L)], mapply(log_bvn_tail, c(-8, -9), c(-7, 2), c(-.3, .8)),
    tolerance = 4e-11)
})

test_that('tail probabilities preserve independent first-score references', {
  p = as.data.frame(bvn_tail_fixture)
  p = p[abs(p$r) <= .95 & p$region != 'threshold_above', ]
  truth = as.matrix(p[c('dlog_a', 'dlog_b', 'dlog_r')])
  score = probit_re_probit_re_rows(p$a, p$b, rep(1, nrow(p)), p$r, 1L)$gradient
  expect_lt(max(abs(score-truth)/pmax(1, abs(truth))), 5e-10)
})

test_that('angular tails handle opposite thresholds and domain boundaries', {
  p = expand.grid(a = c(-8, -30), shift = c(-1e-8, 0, 1e-8),
    r = c(-.95000001, -.95, -.9, .5, .95, .95000001))
  p$b = -p$a+p$shift
  truth = mapply(log_bvn_tail, p$a, p$b, p$r)
  expect_lt(max(abs(log_bvn_value(p$a, p$b, p$r)-truth)), 4e-11)
  opposite = log_bvn_tail_angular(-8, 8, -.9)
  expect_true(is.finite(opposite))
  expect_equal(opposite, log_bvn_tail(-8, 8, -.9), tolerance = 2e-11)
  # Outside the bounded angular/marginal domain, retain the conditional route.
  a = c(-40.00001, -8, -8); b = c(-7, -7, 8.1); r = c(.5, -.95000001, -.5)
  expect_identical(log_bvn_tail_batch(a, b, r), log_bvn_tail_conditional(a, b, r))
})

test_that('the marginal shortcut meets its log-error bound around the switch', {
  a = rep(-8, 3L); r = rep(.5, 3L)
  target = c(.99e-14, 1e-14, 1.01e-14)
  b = r*a+sqrt(1-r*r)*qnorm(-target, log.p = TRUE)
  actual_bound = -pnorm((b-r*a)/sqrt(1-r*r), log.p = TRUE)
  value = log_bvn_value(a, b, r)
  truth = mapply(log_bvn_tail, a, b, r)
  expect_lt(max(abs(value-truth)), 4e-11)
  use = actual_bound <= 1e-14
  expect_true(any(use) && any(!use))
  expect_identical(value[use], pnorm(a[use], log.p = TRUE))
})

test_that('refined angular panels accept checked subset and full-sample tails', {
  # Independent conditional MPFR integrals, 256/384 bits and reversed argument
  # order; first scores use exact partials divided by that reference probability.
  # These literals use the independently evaluated decimal coordinates.
  # Expectations are not production output.
  p = data.frame(
    a = c(-7.0921559672594601, 0.039422852584760401,
          -10.793602239633501, -8.9017929118200492),
    b = c(1.1600070649057801, -8.1707990184018193,
          7.6397448512488504, -2.8343977135260299),
    r = c(-0.162311122382887, 0.162311122382887,
          -0.70188444765944702, 0.70188444765944702),
    log_p = c(-28.750103759643228, -36.498607574776543,
              -62.250273302687098, -42.738394521442132))
  truth = matrix(c(
    7.3582901003717707, 0.81704571713908969, 5.9073492980783078,
    0.16419155942988686, 8.2634752066875006, 1.3984882728673984,
    11.661455519795418, 1.1156882776741019, 12.1425848125574,
    9.011455867356835, 3.7744463745375162e-06, 5.1959003477241087e-05),
    ncol = 3L, byrow = TRUE)
  # Direct acceptance matters: a correct conditional fallback would conceal
  # regression to the coarse angular panels that declined these cases.
  angular = log_bvn_tail_angular(p$a, p$b, p$r)
  expect_true(all(is.finite(angular)))
  expect_lt(max(abs(angular-p$log_p)), 2e-11)
  rows = probit_re_probit_re_rows(p$a, p$b, rep(1, nrow(p)), p$r, 1L)
  expect_lt(max(abs(rows$value-p$log_p)), 2e-11)
  expect_lt(max(abs(rows$gradient-truth)/pmax(1, abs(truth))), 2e-10)
  zero = probit_re_probit_re_rows(p$a, -p$b, rep(0, nrow(p)), -p$r, 1L)
  zero_score = sweep(zero$gradient, 2L, c(1, -1, -1), '*')
  expect_lt(max(abs(zero$value-p$log_p)), 2e-11)
  expect_lt(max(abs(zero_score-truth)/pmax(1, abs(truth))), 2e-10)
})

test_that('centered angular exponents preserve zero limits and mode transitions', {
  # Independent conditional-normal references: 256/384-bit Rmpfr integrals,
  # different GL orders, checked against the swapped 384-bit integral.
  # No expectation comes from PanelSelect.
  # A=0, nearby A>0, then both sides of the constrained-mode transition.
  p = data.frame(
    a = c(-8, -8, -12, -12),
    b = c(8, 7.9999999999990905, 8.3999999998999986, 8.4000000000999986),
    r = c(-.95, -.95, -.7, -.7),
    log_p = c(-35.228519185746308, -35.228519185747224,
      -76.169849168794556, -76.169849168557363))
  score = rbind(
    c(9.0622502630957076, 1.0079505932258048, 5.6637190324786157),
    c(9.0622502630988642, 1.007950593229169, 5.6637190325049405),
    c(12.906918106185822, 1.1859643750766866, 14.420407329804307),
    c(12.906918106008474, 1.1859643748218212, 14.42040732638389))

  angular = log_bvn_tail_angular(p$a, p$b, p$r)
  expect_true(all(is.finite(angular)))
  expect_lt(max(abs(angular-p$log_p)), 2e-11)
  one = probit_re_probit_re_rows(p$a, p$b, rep(1, 4), p$r, 1L)
  expect_lt(max(abs(one$value-p$log_p)), 2e-11)
  expect_lt(max(abs(one$gradient-score)/pmax(1, abs(score))), 2e-10)
  zero = probit_re_probit_re_rows(p$a, -p$b, rep(0, 4), -p$r, 1L)
  zero_score = sweep(zero$gradient, 2L, c(1, -1, -1), '*')
  expect_lt(max(abs(zero$value-p$log_p)), 2e-11)
  expect_lt(max(abs(zero_score-score)/pmax(1, abs(score))), 2e-10)
})
