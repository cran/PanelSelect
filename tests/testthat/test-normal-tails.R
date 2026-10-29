# 1024-bit references for log Phi and its derivatives; see fixtures/README.md.
normal_tail_fixture = dget(testthat::test_path('fixtures', 'normal-tail.R'))$points

normal_tail_error = function(actual, expected) max(abs(actual - expected)/abs(expected))

test_that('log normal derivatives agree with high-precision references in every region', {
  p = normal_tail_fixture
  d = normal_log_derivatives(p$x, 3L)
  deep = p$x < -4
  expect_lt(normal_tail_error(d$value, p$value), 1e-15)
  # The first derivative keeps its log-scale ratio, whose subtraction of large
  # logarithms limits it to about 3e-14 between -30 and -4.
  expect_lt(normal_tail_error(d$first, p$first), 1e-13)
  # Below -4 the continued fraction removes the cancellation (formerly up to
  # 1e-5 near -30). Above -4 the closed forms lose at most a few digits near -4.
  expect_lt(normal_tail_error(d$second[deep], p$second[deep]), 1e-14)
  expect_lt(normal_tail_error(d$third[deep], p$third[deep]), 1e-14)
  expect_lt(normal_tail_error(d$second[!deep], p$second[!deep]), 1e-13)
  expect_lt(normal_tail_error(d$third[!deep], p$third[!deep]), 1e-11)
})

test_that('inner Mills ratios and Taylor derivatives share the tail evaluation', {
  p = normal_tail_fixture
  deep = p$x < -4
  inner = pln_inner_mills(p$x)
  expect_lt(normal_tail_error(inner$ratio, p$mills), 1e-13)
  expect_lt(normal_tail_error(inner$curvature[deep], -p$second[deep]), 1e-14)
  expect_identical(normal_mills(p$x), normal_log_derivatives(p$x, 1L)$first)
  taylor = t(vapply(p$x, inner_taylor_mills_derivatives, numeric(5)))
  reference = cbind(p$first, p$second, p$third, p$fourth, p$fifth)
  expect_lt(normal_tail_error(taylor[deep, 2:5], reference[deep, 2:5]), 1e-14)
  # Above -4 the scalar recurrence is unchanged (largest error at -4).
  expect_lt(normal_tail_error(taylor[!deep, 2:3], reference[!deep, 2:3]), 1e-11)
  expect_lt(normal_tail_error(taylor[!deep, 4:5], reference[!deep, 4:5]), 1e-9)
})

test_that('tail derivatives are continuous across the continued-fraction switch', {
  p = normal_tail_fixture
  take = which(abs(p$x + 4) <= 2e-9 & p$x != -4)
  expect_length(take, 2L)
  d = normal_log_derivatives(p$x[take], 3L)
  # The jump across -4 equals the true change, up to closed-form roundoff.
  expect_lt(abs(diff(d$second) - diff(p$second[take]))/abs(p$second[take[2L]]), 1e-13)
  expect_lt(abs(diff(d$third) - diff(p$third[take]))/abs(p$third[take[2L]]), 1e-11)
})
