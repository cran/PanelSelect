history_guard_fit = function() {
  id = rep(1:8, each = 3L)
  data = data.frame(id = id, x = sin(seq_along(id)),
    y = c(0, 1, 2, 1, 0, 3, 2, 2, 1, 0, 0, 1, 4, 2, 3, 1, 1, 0, 2, 0, 1, 3, 1, 2))
  list(data = data, fit = panel_select(y ~ x, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'gh', nodes = 5L, maxit = 0L)))
}

test_that('a supplied history is never ignored silently', {
  case = history_guard_fit()
  expect_error(predict(case$fit, case$data, history = case$data),
    "conditioning = 'history'")
  expect_error(predict(case$fit, case$data, type = 'link', conditioning = 'history',
    history = case$data), 'Link predictions do not depend on history')
  # The default calls are unaffected.
  expect_length(fitted(case$fit), nrow(case$data))
  expect_length(residuals(case$fit), nrow(case$data))
  # An empty history still means no conditioning information.
  expect_identical(predict(case$fit, case$data, conditioning = 'history',
    history = case$data[0, ]), predict(case$fit, case$data))
})

test_that("history prediction preserves numeric, factor, and character identities", {
  data = data.frame(id = rep(1:3, each = 2), x = c(-1, 0, .2, .5, 1, -.4),
    d = c(1, 0, 1, 1, 0, 1), y = c(2, NA, 0, 3, NA, 1))
  fixed = c("selection.(Intercept)" = .2, selection.x = .1,
    "outcome.(Intercept)" = .3, outcome.x = .2,
    sd_re_selection = .4, sd_re_outcome = .6, cor_re = .25)
  control = panel_control(quadrature = "gh", nodes = 5)
  fit = panel_select(y ~ x, d ~ x, data, "id", "poisson", fixed = fixed, control = control)
  expected = predict(fit, conditioning = "history", history = data, se.fit = TRUE)
  ids = list(1e15 + data$id, factor(data$id, levels = 3:1), paste0("person-", data$id))
  for (id in ids) {
    current = data; current$id = id
    fit = panel_select(y ~ x, d ~ x, current, "id", "poisson", fixed = fixed, control = control)
    expect_equal(predict(fit, conditioning = "history", history = current, se.fit = TRUE),
      expected, tolerance = 1e-12)
  }
})
