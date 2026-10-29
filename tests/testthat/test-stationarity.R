test_that("stationarity uses the analytic Poisson quadratic at a nonoptimal point", {
  data = data.frame(id = 1:5, y = 0:4)
  formula = stats::as.formula("y ~ 1", env = baseenv())
  fit = panel_select(formula, NULL, data, "id", "poisson",
    start = c("outcome.(Intercept)" = log(1.5)), fixed = c(sd_re_outcome = 0),
    control = panel_control(maxit = 0), vcov = 'observed')
  saved = serialize(fit, NULL)
  value = check_stationarity(fit)
  score = sum(data$y)-nrow(data)*1.5
  information = nrow(data)*1.5
  expect_identical(value$status, "available")
  expect_identical(value$optimizer_converged, fit$optimizer$converged)
  expect_equal(unname(value$score), score, tolerance = 1e-12)
  expect_equal(value$information_eigenvalues, information, tolerance = 1e-6)
  expect_equal(unname(value$newton_step), score/information, tolerance = 1e-6)
  expect_equal(value$max_standardized_step, abs(score)/sqrt(information), tolerance = 1e-6)
  expect_identical(serialize(fit, NULL), saved)

  no_information = fit
  no_information$optimizer$information_computed = FALSE
  expect_match(check_stationarity(no_information)$reason, "not computed")
  expect_null(check_stationarity(no_information)$newton_step)
  indefinite = fit
  indefinite$optimizer$information[1, 1] = -1
  expect_match(check_stationarity(indefinite)$reason, "not positive definite")
  expect_null(check_stationarity(indefinite)$standardized_step)
  unnamed = fit
  dimnames(unnamed$optimizer$information) = NULL
  expect_match(check_stationarity(unnamed)$reason, "named observed information")
})

test_that("stationarity reports all-fixed fits without inventing a stopping test", {
  data = data.frame(id = 1:3, y = c(1, 2, 0))
  formula = stats::as.formula("y ~ 1", env = baseenv())
  fit = panel_select(formula, NULL, data, "id", "poisson",
    fixed = c("outcome.(Intercept)" = 0, sd_re_outcome = 0))
  value = check_stationarity(fit)
  expect_identical(value$status, "all_fixed")
  expect_null(value$newton_step)
  expect_true(is.na(value$max_abs_score))
  expect_error(check_stationarity(list()), "panelselect")
})
