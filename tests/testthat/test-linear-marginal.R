# 512-bit references from the definition r ~ N(0, sigma^2 I + lambda^2 11');
# see fixtures/README.md for the independent reference method.
linear_marginal_fixture = fixture_names(dget(testthat::test_path('fixtures', 'linear-marginal.R')))

linear_marginal_error = function(actual, expected) {
  # Elementwise relative error, with a floor for exact zeros in a column.
  scale = pmax(abs(expected), 1e-300)
  max(abs(actual - expected)/scale)
}

test_that('the Linear outcome marginal agrees with 512-bit references', {
  for (case in linear_marginal_fixture$cases) {
    prepared = prepare_data(y ~ x, d ~ w, as.data.frame(case$data), 'id', 'gaussian')
    out = linear_outcome_marginal(case$theta, prepared, conditional = TRUE)
    colnames(out$gradient) = names(case$theta)
    expect_lt(linear_marginal_error(out$value, case$value), 1e-14)
    expect_lt(linear_marginal_error(out$gradient[, colnames(case$gradient)], case$gradient), 1e-14)
    for (field in c('mean', 'variance', 'scale'))
      expect_lt(linear_marginal_error(out[[field]], case[[field]]), 1e-14)
  }
})
