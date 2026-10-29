covariance_object = function(value = diag(2)) {
  parameters = c("outcome.(Intercept)", "outcome.x")
  dimnames(value) = list(parameters, parameters)
  structure(list(coefficients = setNames(c(0, 0), parameters),
    covariance = list(observed = value), control = list(vcov = "observed"),
    optimizer = list(information_computed = TRUE)), class = "panelselect")
}

test_that("explicit prediction covariance rejects invalid matrices", {
  fit = covariance_object()
  named = function(value) {
    dimnames(value) = dimnames(fit$covariance$observed)
    value
  }
  expect_error(prediction_covariance(fit, named(matrix(c(1, 0, 0.8, 1), 2))), "symmetric")
  # Positive variance for one target does not validate an indefinite matrix.
  indefinite = named(matrix(c(1, 2, 2, 1), 2))
  expect_equal(prediction_se(matrix(c(1, 0), 1), indefinite), 1)
  expect_error(prediction_covariance(fit, indefinite), "positive semidefinite")
  expect_error(prediction_covariance(fit, named(diag(c(1, -1e-6))*1e-18)), "positive semidefinite")
  for (value in list(matrix(NA_real_, 2, 2), matrix(Inf, 2, 2), matrix("a", 2, 2),
                    matrix(1+1i, 2, 2), matrix(TRUE, 2, 2)))
    expect_error(prediction_covariance(fit, named(value)), "finite numeric")
  expect_error(prediction_covariance(fit, diag(2)), "named natural-scale")
})

test_that("singular and roundoff covariance inputs remain supported unchanged", {
  fit = covariance_object()
  for (value in list(diag(c(1, 0)), matrix(0, 2, 2), diag(c(1, -1e-14)),
                    diag(c(1, -1e-14))*1e-18, matrix(c(1, 0, 1e-14, 0), 2))) {
    dimnames(value) = dimnames(fit$covariance$observed)
    before = value
    expect_identical(prediction_covariance(fit, value), value)
    expect_identical(value, before)
    expect_equal(prediction_se(matrix(c(0, 1), 1), value), 0)
  }
})

test_that("internal unavailable covariance still propagates only where needed", {
  fit = covariance_object(matrix(c(1, NA, NA, NA), 2))
  value = prediction_covariance(fit)
  expect_identical(value, fit$covariance$observed)
  expect_equal(prediction_se(matrix(c(1, 0), 1), value), 1)
  expect_true(is.na(prediction_se(matrix(c(0, 1), 1), value)))
  expect_error(prediction_covariance(fit, value), "finite numeric")
})
