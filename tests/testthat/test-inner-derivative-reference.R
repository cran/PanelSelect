test_that('the Taylor-jet product equals explicit pairwise sums', {
  setup = inner_taylor_setup
  explicit = function(a, b) vapply(setup$pairs, function(p) sum(a[p[, 1L]]*b[p[, 2L]]), numeric(1))
  magnitude = function(a, b) vapply(setup$pairs, function(p) sum(abs(a[p[, 1L]]*b[p[, 2L]])), numeric(1))
  for (k in 1:20) {
    a = sin(k*seq_len(setup$n)); b = cos(3*k + seq_len(setup$n))^3
    # Summation order differs from sum(); bound by the standard rounding error.
    expect_true(all(abs(inner_taylor_multiply(a, b) - explicit(a, b)) <=
      8*.Machine$double.eps*magnitude(a, b)))
  }
  a[3L] = Inf; b[5L] = NaN
  expect_identical(inner_taylor_multiply(a, b), explicit(a, b))
})

test_that('production adaptive inner PLN agrees with the reference engine', {
  A = c(-.5, .2, .7); B = c(-.3, 1, -1); y = c(0, 4, 9); gamma = .8; tau = -.4
  for (H in c(7L, 13L)) for (order in 0:3) {
    rule = normal_rule(H, 1L); adaptive = rule; adaptive$quadrature = 'agh'
    reference = local_pln_adaptive(B, y, gamma, order, lgamma(y + 1), rule, A, tau, mode_tol = 1e-12)
    actual = integrate_pln_inner(A, B, y, gamma, tau, order, lgamma(y + 1), adaptive,
      probit_re_pln_re_inner_kernel)
    expect_lt(max(abs(actual$value - reference$value)), 1e-14)
    if (order) expect_lt(max(abs(actual$gradient -
      reference$gradient[, seq_len(ncol(actual$gradient)), drop = FALSE])), 1e-13)
  }
  rule = normal_rule(9L, 1L); adaptive = rule; adaptive$quadrature = 'agh'
  reference = local_pln_adaptive(B, y, gamma, 3L, lgamma(y + 1), rule, mode_tol = 1e-12)
  actual = integrate_pln_inner(NULL, B, y, gamma, 0, 3L, lgamma(y + 1), adaptive, pln_re_inner_kernel)
  expect_lt(max(abs(actual$value - reference$value)), 1e-14)
  expect_lt(max(abs(actual$gradient - reference$gradient)), 1e-13)
  expect_lt(max(abs(actual$third - reference$third)), 1e-12)
})
