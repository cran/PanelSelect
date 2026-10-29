test_that('grouped reductions preserve empty panels and weighted design axes', {
  local_test_seed(614)
  ptr = c(0, 0, 1, 4, 4, 7)
  for (p in c(0L, 1L, 4L)) {
    x = matrix(rnorm(7*p), 7, p)
    weights = matrix(rnorm(21), 7, 3)
    saved = serialize(list(x, weights, ptr), NULL)
    expect_equal(group_sum(x, ptr), group_sum_r(x, ptr), tolerance = 1e-13)
    expect_equal(group_weighted_design(x, weights, ptr),
      group_weighted_design_r(x, weights, ptr), tolerance = 1e-13)
    expect_identical(serialize(list(x, weights, ptr), NULL), saved)
  }
  expect_identical(dim(group_sum(matrix(numeric(), 0, 2), 0)), c(0L, 2L))
  expect_identical(dim(group_weighted_design(matrix(numeric(), 0, 2),
    matrix(numeric(), 0, 3), c(0, 0))), c(3L, 3L))
  x = matrix(c(Inf, 1, NA_real_), 3, 1)
  weights = matrix(c(0, 2, 1), 3, 1)
  expect_equal(group_weighted_design(x, weights, c(0, 1, 3)),
    group_weighted_design_r(x, weights, c(0, 1, 3)))
})

test_that('panel mixture weights integrate scores before expanding the design', {
  local_test_seed(731)
  ptr = c(0, 2, 5)
  x = matrix(rnorm(15), 5, 3)
  score = matrix(rnorm(20), 5, 4)
  terms = rnorm(8)
  rule = logsum_weights(terms, 2L)
  expect_equal(rule, logsum_weights_r(terms, 2L), tolerance = 1e-13)
  expect_equal(group_weighted_design(x, score, ptr, rule$weights),
    group_weighted_design_r(x, score, ptr, rule$weights), tolerance = 1e-13)
  score[1:2, 1] = Inf
  rule$weights[1, 1] = 0
  expect_equal(group_weighted_design(x, score, ptr, rule$weights),
    group_weighted_design_r(x, score, ptr, rule$weights), tolerance = 1e-13)
  logs = c(-Inf, NaN, Inf, -1000, -Inf, 0, 0, 0)
  expect_equal(logsum_weights(logs, 4L), logsum_weights_r(logs, 4L))
})

test_that('group boundaries are validated before indexing', {
  x = matrix(1, 3, 2)
  for (ptr in list(numeric(), c(1, 3), c(0, 2), c(0, 4, 3),
                   c(0, 2, 1, 3), c(0, 1.5, 3), c(0, NA, 3), c(0, Inf, 3)))
    expect_error(group_sum(x, ptr), 'boundar')
  expect_error(group_weighted_design(x, matrix(1, 2, 1), c(0, 3)), 'matching rows')
})

test_that('centered mixture moments match R reductions through order three', {
  local_test_seed(915)
  for (n in c(1L, 4L)) for (q in 0:2) for (order in 1:3) {
    d = 4L; nodes = 5L; m = n*nodes
    value = rnorm(m)*12
    g = matrix(rnorm(m*d), m, d)
    h = if (order >= 2L) array(rnorm(m*q*d), c(m, q, d)) else NULL
    t = if (order >= 3L) array(rnorm(m*q*q*d), c(m, q, q, d)) else NULL
    saved = serialize(list(value, g, h, t), NULL)
    expect_equal(mixture_moments(value, g, h, t, n, q),
      mixture_moments_r(value, g, h, t, n, q), tolerance = 1e-12)
    expect_identical(serialize(list(value, g, h, t), NULL), saved)
  }
})

test_that('compact mixture derivatives agree with differentiated R log sums', {
  make_terms = function(theta) {
    variables = lapply(seq_along(theta), function(k) {
      x = jet(theta[k], 4L, 2L, 3L)
      x$gradient[1L, k] = 1
      x
    })
    node = rep(c(-1.4, -.2, .8, 1.7), each = 2L)
    row = rep(c(.7, 1.2), 4L)
    eta = variables[[1L]] + variables[[2L]]*row + .2*tanh(variables[[4L]])*node
    -.5*((node-eta)/exp(variables[[3L]]))^2 - variables[[3L]] +
      .03*variables[[1L]]*variables[[2L]]*variables[[3L]]*node
  }
  theta = c(.2, -.1, .3, -.25)
  terms = make_terms(theta)
  observed = mixture_moments(terms$value, terms$gradient, terms$hessian, terms$third, 2L, 2L)
  reference = jet_log_mixture_r(terms, 2L)
  for (field in names(observed)) expect_equal(observed[[field]], reference[[field]], tolerance = 1e-11)
  for (j in 1:2) {
    step = rep(0, 4L); step[j] = 1e-5
    above = jet_log_mixture_r(make_terms(theta+step), 2L)
    below = jet_log_mixture_r(make_terms(theta-step), 2L)
    expect_equal(matrix(observed$third[, j, , ], 2L, 8L),
      matrix((above$hessian-below$hessian)/(2e-5), 2L, 8L), tolerance = 2e-8)
  }
  shifted = mixture_moments(terms$value+rep(c(-1000, 1000), 4L), terms$gradient,
    terms$hessian, terms$third, 2L, 2L)
  expect_equal(shifted$value, observed$value+c(-1000, 1000), tolerance = 1e-12)
  expect_equal(shifted$third, observed$third, tolerance = 1e-10)
})

test_that('mixture moments distinguish invalid terms, absent mass and underflow', {
  value = c(-Inf, NaN, Inf, -1000, -Inf, 0, 0, 0)
  g = matrix(c(Inf, 1, 1, Inf, Inf, 2, 2, 3), 8, 1)
  h = array(g, c(8, 1, 1)); t = array(g, c(8, 1, 1, 1))
  out = mixture_moments(value, g, h, t, 4L, 1L)
  expect_equal(out, mixture_moments_r(value, g, h, t, 4L, 1L))
  expect_equal(out$value[c(1, 4)], c(-Inf, 0))
  expect_equal(out$gradient[c(1, 4), 1], c(0, 3))
  expect_true(all(is.nan(out$third[c(2, 3), , , ])))
  expect_equal(logsum_score(value, g, 4L), logsum_score_r(value, g, 4L))
  empty = mixture_moments(numeric(), matrix(numeric(), 0, 2),
    array(numeric(), c(0, 1, 2)), array(numeric(), c(0, 1, 1, 2)), 0, 1)
  expect_identical(dim(empty$third), c(0L, 1L, 1L, 2L))
  absent = mixture_moments(numeric(), matrix(numeric(), 0, 0), n = 3, q = 0)
  expect_identical(absent$value, rep(-Inf, 3L))
})

test_that('empty numerical blocks keep shapes without constructing null subviews', {
  ptr = c(0, 2, 5)
  expect_identical(group_sum(matrix(numeric(), 5, 0), ptr), matrix(numeric(), 2, 0))
  weights = matrix(seq_len(10), 5, 2)
  expect_equal(group_weighted_design(matrix(numeric(), 5, 0), weights, ptr),
    group_weighted_design_r(matrix(numeric(), 5, 0), weights, ptr))
  empty = mixture_moments(numeric(), matrix(numeric(), 0, 2),
    array(numeric(), c(0, 1, 2)), array(numeric(), c(0, 1, 1, 2)), n = 3, q = 1)
  expect_equal(empty, mixture_moments_r(numeric(), matrix(numeric(), 0, 2),
    array(numeric(), c(0, 1, 2)), array(numeric(), c(0, 1, 1, 2)), n = 3, q = 1))
  x = matrix(c(1, Inf, 2, 3, 4), 5, 1)
  expect_equal(group_weighted_design(x, matrix(numeric(), 5, 0), ptr, matrix(numeric(), 2, 0)),
    group_weighted_design_r(x, matrix(numeric(), 5, 0), ptr, matrix(numeric(), 2, 0)))
})

test_that('mixture derivative dimensions cannot be silently recycled', {
  g = matrix(1, 6, 3)
  expect_error(mixture_moments(1:6, g, n = 2.5), 'nonnegative integer')
  expect_error(mixture_moments(1:6, g, n = 2, q = 4), 'exceeds')
  expect_error(mixture_moments(1:6, g, array(1, c(6, 2, 2)), n = 2, q = 2), 'dimensions')
  expect_error(mixture_moments(1:6, g, third = array(1, c(6, 2, 2, 3)), n = 2, q = 2), 'require')
})

test_that('batched affine maps preserve node order and analytic derivatives', {
  local_test_seed(615)
  for (groups in c(1L, 3L)) for (q in 0:2) for (p in c(0L, 1L, 4L)) {
    base = matrix(rnorm(5*q), 5, q)
    location = matrix(rnorm(groups*q), groups, q)
    transform = array(rnorm(groups*q*q), c(groups, q, q))
    dm = array(rnorm(groups*q*p), c(groups, q, p))
    dl = array(rnorm(groups*q*q*p), c(groups, q, q, p))
    saved = serialize(list(base, location, transform, dm, dl), NULL)
    out = affine_batch(base, location, transform, dm, dl)
    expect_equal(out, affine_batch_r(base, location, transform, dm, dl), tolerance = 1e-13)
    expect_identical(serialize(list(base, location, transform, dm, dl), NULL), saved)
    score = matrix(rnorm(groups*5*q), groups*5, q)
    expect_equal(row_contract(score, out$gradient), row_contract_r(score, out$gradient), tolerance = 1e-13)
    plain = affine_batch(base, location, transform)
    expect_identical(names(plain), c('nodes', 'gradient'))
    expect_null(plain$gradient)
    expect_equal(plain$nodes, out$nodes)
    if (p && q) {
      shift_m = matrix(dm[, , 1L], groups, q)
      shift_l = array(dl[, , , 1L], c(groups, q, q))
      above = affine_batch(base, location+1e-5*shift_m, transform+1e-5*shift_l)
      below = affine_batch(base, location-1e-5*shift_m, transform-1e-5*shift_l)
      expect_equal(matrix(out$gradient[, , 1L], groups*5, q),
        (above$nodes-below$nodes)/(2e-5), tolerance = 1e-10)
    }
  }
})

test_that('affine and contraction dimensions are validated', {
  base = matrix(1, 3, 2); location = matrix(1, 4, 2)
  transform = array(1, c(4, 2, 2))
  expect_error(affine_batch(base, location, array(1, c(4, 2, 3))), 'dimensions')
  expect_error(affine_batch(base, location, transform, array(1, c(4, 2, 1))), 'Both')
  expect_error(row_contract(matrix(1, 3, 2), array(1, c(3, 3, 1))), 'dimensions')
  score = matrix(c(Inf, 1, 2, 3), 2, 2)
  derivative = array(0, c(2, 2, 1))
  expect_equal(row_contract(score, derivative), row_contract_r(score, derivative))
})
