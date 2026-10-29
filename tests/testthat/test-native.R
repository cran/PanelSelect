native_mixture_reference = function(values, score, n) {
  terms = matrix(values, n)
  answer = matrix(0, n, ncol(score)+1L)
  for (i in seq_len(n)) {
    index = i+n*(seq_len(ncol(terms))-1L)
    x = terms[i, ]
    if (anyNA(x) || any(x == Inf)) { answer[i, ] = NaN; next }
    if (all(x == -Inf)) { answer[i, 1] = -Inf; next }
    high = max(x); mass = exp(x-high)
    answer[i, 1] = high+log(sum(mass))
    keep = mass > 0
    answer[i, -1L] = drop(crossprod(mass[keep]/sum(mass), score[index[keep], , drop = FALSE]))
  }
  answer
}

test_that('native mixture values and scores agree without mutating inputs', {
  local_test_seed(81L)
  for (n in c(1L, 3L, 17L)) for (k in c(1L, 5L)) for (p in c(0L, 1L, 6L)) {
    values = rnorm(n*k)*1000
    score = matrix(rnorm(n*k*p), n*k, p)
    before = serialize(list(values, score), NULL)
    expect_equal(logsum_score_cpp(values, score, n),
      native_mixture_reference(values, score, n), tolerance = 1e-12)
    expect_identical(serialize(list(values, score), NULL), before)
  }
})

test_that('native mixture distinguishes absent mass and invalid components', {
  values = c(-Inf, NaN, Inf, -Inf, 0.3, 1000, 0, -Inf)
  score = matrix(c(Inf, 1, 1, Inf, 2, 1, 1, Inf), 8, 1)
  expect_equal(logsum_score_cpp(values, score, 4L), native_mixture_reference(values, score, 4L))
  expect_identical(dim(logsum_score_cpp(numeric(), matrix(numeric(), 0, 3), 0L)), c(0L, 4L))
  expect_error(logsum_score_cpp(1:3, matrix(0, 3, 1), 2L), 'must agree')
  expect_error(logsum_score_cpp(1:3, matrix(0, 2, 1), 3L), 'must agree')
  for (rows in c(-1, 1.5, NA_real_, Inf, .Machine$integer.max+1))
    expect_error(logsum_score_cpp(1:3, matrix(0, 3, 1), rows), 'nonnegative integer')
})
