derivative_example = function(values, order = 3L) {
  x = seed_jets(values, order = order)
  a = x[[1]]*exp(x[[3]]) + x[[2]]^2
  b = x[[2]]/(2+x[[1]]^2) + x[[3]]
  r = tanh(x[[3]])*0.6
  jet_log_pnorm(a) + jet_log_bvn(a, b, r) + exp(a/5)*sqrt(b*b+2) +
    tanh(a*b/4) + jet_logadd(-a*a, -b*b)
}

test_that('mixed derivative composition agrees with perturbed evaluations', {
  values = c(0.3, -0.4, 0.25)
  actual = derivative_example(values)
  h = 1e-5
  for (j in seq_along(values)) {
    plus = minus = values; plus[j] = plus[j]+h; minus[j] = minus[j]-h
    upper = derivative_example(plus); lower = derivative_example(minus)
    expect_equal(actual$gradient[1, j], (upper$value-lower$value)/(2*h), tolerance = 1e-7)
    if (j <= 2L) {
      expect_equal(drop(actual$hessian[1, j, ]),
        drop((upper$gradient-lower$gradient)/(2*h)), tolerance = 1e-7)
    }
    for (k in 1:2) for (l in 1:2) {
      expect_equal(actual$third[1, k, l, j],
        (upper$hessian[1, k, l]-lower$hessian[1, k, l])/(2*h), tolerance = 2e-6)
    }
  }
  for (order in 0:2) {
    reduced = derivative_example(values, order)
    expect_equal(reduced$value, actual$value)
    if (order >= 1) expect_equal(reduced$gradient, actual$gradient)
    if (order >= 2) expect_equal(reduced$hessian, actual$hessian)
  }
})

test_that('bivariate CDF retains correlation scores at zero', {
  x = seed_jets(c(0.3, -0.4, 0))
  out = jet_log_bvn(x[[1]], x[[2]], x[[3]])
  mills = dnorm(c(0.3, -0.4))/pnorm(c(0.3, -0.4))
  expect_equal(out$value, sum(pnorm(c(0.3, -0.4), log.p = TRUE)), tolerance = 1e-13)
  expect_equal(out$gradient[1, ], c(mills, prod(mills)), tolerance = 1e-12)
  expect_gt(abs(out$gradient[1, 3]), 0.1)
})

test_that('bivariate tail derivatives remain finite and match differences', {
  for (point in list(c(-12, -11, 0.5), c(-7, -6, -0.4), c(0.2, -0.3, 0.9))) {
    x = seed_jets(point)
    actual = jet_log_bvn(x[[1]], x[[2]], x[[3]])
    expect_true(all(is.finite(unlist(actual[c('value', 'gradient', 'hessian', 'third')]))))
    h = 1e-5
    for (j in seq_along(point)) {
      plus = minus = point; plus[j] = plus[j]+h; minus[j] = minus[j]-h
      expected = (log_bvn_value(plus[1], plus[2], plus[3])-
        log_bvn_value(minus[1], minus[2], minus[3]))/(2*h)
      expect_equal(actual$gradient[1, j], expected, tolerance = 3e-5)
      p = seed_jets(plus); m = seed_jets(minus)
      upper = jet_log_bvn(p[[1]], p[[2]], p[[3]])
      lower = jet_log_bvn(m[[1]], m[[2]], m[[3]])
      for (k in 1:2) for (l in 1:2)
        expect_equal(actual$third[1, k, l, j],
          (upper$hessian[1, k, l]-lower$hessian[1, k, l])/(2*h), tolerance = 3e-4)
    }
  }
})

test_that('normal-tail derivative expansion is smooth and stable', {
  x = jet(c(-100, -40, -30-1e-6, -30+1e-6, 0, 8), 1L, 1L, 3L)
  x$gradient[, 1] = 1
  out = jet_log_pnorm(x)
  expect_true(all(is.finite(unlist(out[c('value', 'gradient', 'hessian', 'third')]))))
  expect_equal(out$gradient[, 1], exp(dnorm(x$value, log = TRUE)-out$value), tolerance = 1e-10)
  expect_equal(out$hessian[3, 1, 1], out$hessian[4, 1, 1], tolerance = 1e-8)
  expect_lt(abs(out$third[3, 1, 1, 1]-out$third[4, 1, 1, 1]), 1e-9)
  expect_equal(out$hessian[1, 1, 1], -0.99990005995, tolerance = 1e-11)
})

test_that('broadcast, row operations, and empty groups preserve derivatives', {
  parameter = jet(0.4, 1L, 1L, 3L); parameter$gradient[1, 1] = 1
  x = exp(parameter*c(1, 2, 3))
  expect_equal(x$value, exp(0.4*c(1, 2, 3)))
  expect_equal(x$gradient[, 1], c(1, 2, 3)*x$value)
  expect_equal(x$third[, 1, 1, 1], c(1, 2, 3)^3*x$value)
  grouped = jet_group_sum(x, c(1L, 3L, 1L), 3L)
  expect_equal(grouped$value, c(sum(x$value[c(1, 3)]), 0, x$value[2]))
  expect_equal(grouped$third[, 1, 1, 1], c(sum(x$third[c(1, 3), 1, 1, 1]), 0, x$third[2, 1, 1, 1]))
  selected = jet_where(c(TRUE, FALSE, TRUE), x, 0)
  expect_equal(selected$gradient[, 1], x$gradient[, 1]*c(1, 0, 1))
  broadcast = jet_where(c(TRUE, FALSE, TRUE), parameter, 0)
  expect_equal(broadcast$gradient[, 1], c(1, 0, 1))
  expect_equal(jet_subset(x, 2L)$third[1, 1, 1, 1], x$third[2, 1, 1, 1])
  expect_equal(jet_replace(x, integer(), jet_subset(x, integer()))$value, x$value)
  expect_error(x+c(1, 2), 'matching rows')
})

test_that('log mixture handles absent mass and has exact derivatives', {
  parameter = jet(c(-2, 0.3), 1L, 1L, 3L); parameter$gradient[, 1] = 1
  initial = jet(c(-Inf, -Inf), 1L, 1L, 3L)
  expect_equal(jet_logadd(initial, parameter), parameter)
  initial$gradient[] = initial$hessian[] = initial$third[] = Inf
  expect_equal(jet_logadd(initial, parameter), parameter)
  expect_equal(jet_logadd(parameter, initial), parameter)
  absent = jet_logadd(initial, initial)
  expect_equal(absent$value, c(-Inf, -Inf))
  expect_equal(absent$gradient, matrix(0, 2, 1))
  invalid = jet_logadd(parameter, c(NaN, Inf))
  expect_true(is.nan(invalid$value[1]))
  expect_equal(invalid$value[2], Inf)
  value = jet_logadd(parameter, 0)
  p = plogis(parameter$value)
  expect_equal(value$gradient[, 1], p)
  expect_equal(value$hessian[, 1, 1], p*(1-p))
  expect_equal(value$third[, 1, 1, 1], p*(1-p)*(1-2*p))
})

test_that('zero-dimensional derivative blocks and powers at zero are valid', {
  x = jet(c(0, 1), 0L, 0L, 0L)
  expect_equal((x^2)$value, c(0, 1))
  expect_equal(dim(jet_group_sum(x, c(1L, 1L), 2L)$gradient), c(2L, 0L))
  y = jet(0, 1L, 1L, 3L); y$gradient[1, 1] = 1
  expect_equal((y^2)$third[1, 1, 1, 1], 0)
  expect_equal((y^0)$gradient[1, 1], 0)
})
