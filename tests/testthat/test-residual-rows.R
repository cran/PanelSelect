test_that("residuals reject new-row predictions after normal argument matching", {
  data = data.frame(id = rep(1:3, each = 2), x = c(-1, 0, .2, .5, 1, -.4),
    d = c(1, 0, 1, 1, 0, 1), y = c(2, NA, 0, 3, NA, 1))
  fit = panel_select(y ~ x, d ~ x, data, "id", "poisson",
    fixed = c("selection.(Intercept)" = .2, selection.x = .1,
      "outcome.(Intercept)" = .3, outcome.x = .2,
      sd_re_selection = .4, sd_re_outcome = .6, cor_re = .25),
    control = panel_control(quadrature = "gh", nodes = 3))
  target = fit$data[nrow(fit$data):1L, ]
  expected = residuals(fit)
  expect_equal(residuals(fit, newdata = NULL), expected)
  expect_equal(residuals(fit, newd = NULL), expected)
  expect_equal(residuals(fit, NULL), expected)
  expect_equal(residuals(fit, NULL, "joint"), residuals(fit, type = "joint"))
  expect_equal(residuals(fit, typ = "joint"), residuals(fit, type = "joint"))
  for (arguments in list(list(newdata = target), list(newd = target), list(target)))
    expect_error(do.call(residuals, c(list(object = fit), arguments)), "retained fitting rows")
 
  # Conditioning changes predicted means, but responses still use fitting order.
  history = fit$data[fit$data$id == 1L, ]
  predicted = predict(fit, type = "selected_mean", conditioning = "history", history = history)
  response = fit$data$y; response[fit$data$d == 0] = NA_real_
  expect_equal(residuals(fit, conditioning = "history", history = history), response - predicted)
})
