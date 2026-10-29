test_that('parameter transforms require canonical parameter order', {
  info = data.frame(name = c('outcome.x', 'sd_re_outcome', 'cor_re'),
    type = c('coefficient', 'sd', 'correlation'), fixed = FALSE)
  theta = c(outcome.x = .5, sd_re_outcome = 2, cor_re = .3)
  eta = parameter_to_optimizer(theta, info)
  expect_equal(unname(parameter_from_optimizer(eta, theta, info)), unname(theta), tolerance = 1e-15)
  shuffled = theta[c(2L, 1L, 3L)]
  expect_error(parameter_to_optimizer(shuffled, info), 'parameter order')
  expect_error(parameter_from_optimizer(eta, shuffled, info), 'parameter order')
  expect_error(parameter_jacobian(shuffled, info), 'parameter order')
})

test_that('normal Mills ratios keep missing values in place', {
  value = normal_mills(c(NA, -40, NA, 0))
  expect_true(is.na(value[1L]) && is.na(value[3L]))
  expect_equal(value[4L], 2*stats::dnorm(0), tolerance = 1e-15)
  # log-scale density and probability are accurate at -40, so their ratio is too.
  expect_equal(value[2L], exp(stats::dnorm(-40, log = TRUE) - stats::pnorm(-40, log.p = TRUE)),
    tolerance = 1e-13)
})
