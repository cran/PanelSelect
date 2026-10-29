# Independent materialized finite sums. The reference retains node pairs only
# in these small fixtures; the production reducer streams them by cell label.
bilinear_cells_reference = function(a, b, x, y, cross, score_a, score_b,
                                     cross_score, cell_a, cell_b, ncells_a, ncells_b) {
  p = ncol(score_a)
  out = matrix(0, ncells_a*ncells_b, p+1L); out[, 1L] = -Inf
  terms = outer(a, b, '+')+outer(x, cross*y)
  for (v in seq_len(ncells_b)) for (u in seq_len(ncells_a)) {
    i = which(cell_a == u); j = which(cell_b == v)
    if (!length(i) || !length(j)) next
    log_mass = terms[i, j, drop = FALSE]; high = max(log_mass)
    if (high == -Inf) next
    weights = exp(log_mass-high); total = sum(weights)
    row = u+ncells_a*(v-1L)
    out[row, 1L] = high+log(total)
    weights = weights/total
    for (k in seq_len(p)) {
      local = outer(score_a[i, k], score_b[j, k], '+')+
        outer(cross_score[k]*x[i], y[j])
      local[weights == 0] = 0
      out[row, k+1L] = sum(weights*local)
    }
  }
  out
}

test_that('streamed cell likelihoods and scores equal independent finite tensors', {
  local_test_seed(610073L)
  fixtures = list(
    list(a = c(-3, -.7, -1.2, -4, -2), b = c(-.2, -2, -.8, -3),
      x = c(-2, .1, 3, -.7, 1.8), y = c(-1, 2, .3, -.4), cross = -.35,
      cell_a = c(3, 1, 3, 1, 2), cell_b = c(2, 1, 2, 1), ncells_a = 4, ncells_b = 3),
    list(a = c(-1001, -1003, -1002), b = c(-2, -4, -.3, -8),
      x = c(-4, .25, 1.3), y = c(2, -1, 0, .8), cross = .7,
      cell_a = c(1, 2, 1), cell_b = c(2, 1, 1, 2), ncells_a = 2, ncells_b = 2),
    list(a = c(-.3, -2, -5), b = c(-.4, -Inf, -.8),
      x = c(-3, 0, 4), y = c(-2, 1, .7), cross = 0,
      cell_a = c(1, 1, 3), cell_b = c(1, 2, 1), ncells_a = 3, ncells_b = 3))
  for (z in fixtures) {
    scored = unscored = NULL
    for (p in c(0L, 3L)) {
      z$score_a = matrix(rnorm(length(z$a)*p), length(z$a), p)
      z$score_b = matrix(rnorm(length(z$b)*p), length(z$b), p)
      z$cross_score = rnorm(p)
      saved = serialize(z, NULL)
      actual = do.call(bilinear_cells_cpp, z)
      reference = do.call(bilinear_cells_reference, z)
      expect_equal(actual, reference, tolerance = 3e-13)
      expect_identical(dim(actual), c(as.integer(z$ncells_a*z$ncells_b), p+1L))
      expect_identical(serialize(z, NULL), saved)
      absent = reference[, 1L] == -Inf
      expect_true(all(actual[absent, -1L, drop = FALSE] == 0))
      if (p) scored = actual[, 1L] else unscored = actual[, 1L]
    }
    expect_identical(scored, unscored)
  }
})

test_that('cell IDs are one based, unsorted and u-cell fastest', {
  z = list(a = log(c(.1, .2)), b = log(c(.3, .4)), x = c(-1, 1), y = c(.2, -.3),
    cross = 0, score_a = matrix(numeric(), 2, 0), score_b = matrix(numeric(), 2, 0),
    cross_score = numeric(), cell_a = c(2, 1), cell_b = c(2, 1), ncells_a = 3, ncells_b = 2)
  actual = do.call(bilinear_cells_cpp, z)
  expect_equal(drop(actual), log(c(.08, .04, 0, .06, .03, 0)), tolerance = 1e-14)
  expect_equal(actual, do.call(bilinear_cells_reference, z), tolerance = 1e-14)
})

test_that('zero embedded weights and underflow mask otherwise absent derivatives', {
  z = list(a = c(0, -Inf, -1000), b = c(-.4, -Inf, -.7),
    x = c(-1, .3, 2), y = c(-2, 1, .7), cross = .12,
    score_a = matrix(c(.2, NaN, Inf, -.3, Inf, NaN, .7, -Inf, Inf), 3, 3),
    score_b = matrix(c(.4, Inf, .3, -.5, NaN, .6, .8, -Inf, -.2), 3, 3),
    cross_score = c(.3, -.4, .7), cell_a = c(1, 2, 1), cell_b = c(1, 2, 1),
    ncells_a = 3, ncells_b = 3)
  actual = do.call(bilinear_cells_cpp, z)
  expect_equal(actual, do.call(bilinear_cells_reference, z), tolerance = 3e-13)
  expect_true(is.finite(actual[1L, 1L]))
  expect_true(all(is.finite(actual[1L, -1L])))
  expect_true(all(actual[-1L, 1L] == -Inf))
  expect_true(all(actual[-1L, -1L] == 0))
  z$a[] = -Inf
  absent = do.call(bilinear_cells_cpp, z)
  expect_identical(absent[, 1L], rep(-Inf, 9L))
  expect_true(all(absent[, -1L] == 0))
  z$a = z$x = z$cell_a = numeric(); z$score_a = matrix(numeric(), 0, 3)
  expect_identical(do.call(bilinear_cells_cpp, z), absent)
})

test_that('cell reductions reject malformed dimensions, labels and positive-weight scores', {
  z = list(a = 0, b = 0, x = 0, y = 0, cross = 0,
    score_a = matrix(.2, 1, 1), score_b = matrix(.3, 1, 1), cross_score = .1,
    cell_a = 1, cell_b = 1, ncells_a = 1, ncells_b = 1)
  for (field in c('ncells_a', 'ncells_b')) for (invalid in c(-1, .5, Inf, NA_real_, .Machine$integer.max+1)) {
    trial = z; trial[[field]] = invalid
    expect_error(do.call(bilinear_cells_cpp, trial), 'dimensions')
  }
  for (field in c('cell_a', 'cell_b')) for (invalid in c(0, 2, .5, Inf, NA_real_)) {
    trial = z; trial[[field]] = invalid
    expect_error(do.call(bilinear_cells_cpp, trial), 'cell labels')
  }
  for (field in c('x', 'y', 'cross', 'a', 'b')) {
    trial = z; trial[[field]] = NaN
    expect_error(do.call(bilinear_cells_cpp, trial), 'nonfinite')
  }
  trial = z; trial$score_b = matrix(0, 1, 2)
  expect_error(do.call(bilinear_cells_cpp, trial), 'dimensions')
  trial = z; trial$cross_score = numeric()
  expect_error(do.call(bilinear_cells_cpp, trial), 'dimensions')
  trial = z; trial$cell_a = numeric()
  expect_error(do.call(bilinear_cells_cpp, trial), 'dimensions')
  trial = z; trial$score_a[1L, 1L] = Inf
  expect_error(do.call(bilinear_cells_cpp, trial), 'positive-weight score')
  trial = z; trial$x = trial$y = 2; trial$cross_score = Inf
  expect_error(do.call(bilinear_cells_cpp, trial), 'positive-weight score')
  trial = z; trial$ncells_a = trial$ncells_b = 50000
  expect_error(do.call(bilinear_cells_cpp, trial), 'dimensions')
})
