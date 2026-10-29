subset_order_fixture = function(selection = TRUE) {
  data = data.frame(id = c(20, 10, 30, 10, 30, 10, 40, 40, 40),
    x = c(.2, -.4, .3, .1, -.2, .5, .6, -.5, .4),
    os = c(.1, -.2, .3, .2, -.1, -.2, .1, -.1, .2),
    oy = c(-.1, .2, -.2, .1, .3, -.1, .2, -.2, .1),
    d = c(1, 1, 1, 0, 0, 1, 0, 0, 0),
    y = c(1, 0, 1, NA, NA, 1, NA, NA, NA))
  if (!selection) {
    data$d = 1
    data$y = c(1, 0, 1, 0, 1, 1, 0, 1, 0)
  }
  outcome = stats::as.formula('y ~ x + offset(oy)', env = asNamespace('stats'))
  select = if (selection)
    stats::as.formula('d ~ x + offset(os)', env = asNamespace('stats')) else NULL
  prepared = prepare_data(outcome, select, data, 'id', 'probit', 'fail')
  theta = c('outcome.(Intercept)' = .2, outcome.x = -.3, sd_re_outcome = .7)
  if (selection) theta = c('selection.(Intercept)' = -.15, selection.x = .25,
    theta, sd_re_selection = .6, cor_re = .35, cor_error = 0)
  list(prepared = prepared, theta = theta,
    fixed = if (selection) c(cor_error = 0) else numeric())
}

test_that('subsets retain requested panel order and every row-aligned field', {
  f = subset_order_fixture()
  prepared = prepare_model_data(f$prepared)
  prepared$active = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  prepared$integration_cache = list(stale = TRUE)
  original = serialize(prepared, NULL)
  for (panels in list(c(3L, 1L), c(4L, 2L, 1L), 1:4, c(1L, 3L), 2L,
                      integer(), numeric(), NULL)) {
    actual = subset_panels(prepared, panels)
    rows = as.integer(unlist(lapply(panels, function(i) which(prepared$group == i))))
    lengths = vapply(panels, function(i) sum(prepared$group == i), integer(1))
    expect_identical(actual$group, rep(seq_along(panels), lengths))
    expect_equal(actual$group_ptr, c(0L, cumsum(lengths)))
    expect_identical(actual$id_values, prepared$id_values[panels])
    expect_equal(c(actual$n_rows, actual$n_panels, actual$n_units),
      c(length(rows), length(panels), length(panels)))
    for (name in c('x', 'w'))
      expect_identical(actual[[name]], prepared[[name]][rows, , drop = FALSE])
    for (name in c('y', 'd', 'offset_outcome', 'offset_selection', 'selected', 'active'))
      expect_identical(actual[[name]], prepared[[name]][rows])
    expect_identical(actual$original_index, prepared$original_index[rows])
    expect_identical(actual$keep, sort(prepared$original_index[rows]))
    expect_identical(actual$inverse_index, order(prepared$original_index[rows]))
    expect_identical(actual$design, prepared$design)
    expect_null(actual$model_data)
    expect_null(actual$integration_cache)
    blocks = prepare_model_data(actual)$model_data
    expect_equal(blocks$active$count, tabulate(actual$group[actual$active], actual$n_panels))
    expect_equal(blocks$selected$count,
      tabulate(actual$group[actual$active & actual$d == 1], actual$n_panels))
    expect_equal(blocks$unselected$count,
      tabulate(actual$group[actual$active & actual$d == 0], actual$n_panels))
  }
  expect_identical(serialize(prepared, NULL), original)
})

test_that('unsorted subset likelihoods and scores equal individual panels', {
  for (selection in c(FALSE, TRUE)) {
    f = subset_order_fixture(selection)
    for (method in c('gh', 'agh')) {
      control = panel_control(quadrature = method, nodes = 7)
      for (panels in list(c(3L, 1L), c(4L, 2L, 1L), c(2L, 4L), 1:4, 3L)) {
        subset = subset_panels(f$prepared, panels)
        actual = evaluate_model(f$theta, subset, control, f$fixed)
        expected = lapply(panels, function(i)
          evaluate_model(f$theta, subset_panels(f$prepared, i), control, f$fixed))
        expect_equal(unname(actual$loglik), unname(unlist(lapply(expected, `[[`, 'loglik'))),
          tolerance = 1e-11)
        expect_equal(unname(actual$score), unname(do.call(rbind, lapply(expected, `[[`, 'score'))),
          tolerance = 1e-10)
        expect_equal(evaluate_model(f$theta, subset, control, f$fixed, score = FALSE)$loglik,
          actual$loglik, tolerance = 1e-11)
      }
    }
  }
})

test_that('invalid panel indices reject before selection or recycling', {
  f = subset_order_fixture()
  bad = list(c(1, 1), c(1, NA), NA_integer_, NaN, Inf, -Inf, 0, -1, 5,
    1.2, c(1, 2.5), '1', TRUE, list(1), matrix(1, 1, 1), 1+0i)
  for (panels in bad)
    expect_error(subset_panels(f$prepared, panels), 'distinct finite integer indices')
  expect_equal(subset_panels(f$prepared, c(third = 3, first = 1))$id_values, c(30, 10))
})

test_that('zero-row groups and minimal prepared objects remain aligned', {
  f = subset_order_fixture()
  prepared = f$prepared
  prepared$group = c(1L, 1L, 1L, 3L, 4L, 4L, 5L, 5L, 5L)
  prepared$n_panels = prepared$n_units = 5L
  prepared$id_values = c(10, 15, 20, 30, 40)
  prepared$group_ptr = c(0L, 3L, 3L, 4L, 6L, 9L)
  subset = subset_panels(prepared, c(4L, 2L, 1L))
  expect_equal(subset$group, c(1L, 1L, 3L, 3L, 3L))
  expect_equal(subset$group_ptr, c(0L, 2L, 2L, 5L))
  expect_equal(subset$id_values, c(30, 15, 10))
  expect_equal(prepare_model_data(subset)$model_data$active$count, c(2L, 0L, 3L))
  prepared$original_index = prepared$keep = prepared$inverse_index = NULL
  expect_null(subset_panels(prepared, c(4L, 2L))$original_index)
})
