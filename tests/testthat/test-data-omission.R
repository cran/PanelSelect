# Omitted rows must not leave design columns for levels that no retained row
# uses. The expectation is the same model on user-subsetted data, evaluated at
# its automatic starts (maxit = 0): designs, starts and likelihood must agree.
omission_fixture = function() {
  id = rep(1:8, each = 3L)
  data = data.frame(id = id, x = sin(seq_along(id)), x2 = cos(2*seq_along(id)),
    f = factor(rep(c('a', 'b', 'c'), length.out = length(id))),
    y = c(0, 1, 2, 1, 0, 3, 2, 2, 1, 0, 0, 1, 4, 2, 3, 1, 1, 0, 2, 0, 1, 3, 1, 2))
  data$x[data$f == 'c'] = NA
  data
}

test_that('a factor level present only in omitted rows is dropped, as in a pre-subsetted fit', {
  data = omission_fixture()
  complete = droplevels(data[data$f != 'c', ])
  reference = panel_select(y ~ x + f, NULL, complete, 'id', 'poisson',
    control = panel_control(quadrature = 'gh', nodes = 5L, maxit = 0L))
  for (action in c('omit', 'exclude')) {
    fit = panel_select(y ~ x + f, NULL, data, 'id', 'poisson',
      control = panel_control(quadrature = 'gh', nodes = 5L, maxit = 0L), na_action = action)
    expect_identical(coef(fit), coef(reference))
    expect_identical(fit$loglik, reference$loglik)
    expect_identical(fit$prepared$design$outcome$xlevels, list(f = c('a', 'b')))
  }
  excluded = panel_select(y ~ x + f, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'gh', nodes = 5L, maxit = 0L), na_action = 'exclude')
  expect_length(residuals(excluded), nrow(data))
  expect_true(all(is.na(residuals(excluded)[data$f == 'c'])))
  expect_error(predict(excluded, data.frame(id = 1, x = 0, x2 = 0, f = 'c')), 'new level')
})

test_that('data-dependent terms keep their full-data evaluation when levels are dropped', {
  data = omission_fixture()
  fit = panel_select(y ~ poly(x2, 2) + f + x, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'gh', nodes = 5L, maxit = 0L), na_action = 'omit')
  # Bases are evaluated on every supplied row before omission, as before.
  basis = stats::poly(data$x2, 2)[fit$prepared$original_index, , drop = FALSE]
  columns = grep('^poly', colnames(fit$prepared$x))
  expect_equal(unname(fit$prepared$x[, columns]), unname(basis), tolerance = 0)
  expect_false(any(colSums(abs(fit$prepared$x)) == 0))
})

test_that('a linear outcome without selection reports the missing selection equation', {
  data = data.frame(id = rep(1:3, each = 2L), x = 1:6, y = c(.1, .2, .4, .3, .2, .5))
  expect_error(panel_select(y ~ x, NULL, data, 'id', 'gaussian'), 'requires a selection equation')
})
