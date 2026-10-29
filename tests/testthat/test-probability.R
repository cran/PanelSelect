test_that("the bivariate probability respects normalization and known identities", {
  correlations = c(-.95, -.4, 0, .6, .95)
  expect_equal(exp(log_bvn_value(0, 0, correlations)),
               .25 + asin(correlations)/(2*pi), tolerance = 1e-13)
  for (r in correlations) {
    a = -.7; b = 1.2
    p = exp(log_bvn_value(a, b, r))
    expect_equal(p + exp(log_bvn_value(a, -b, -r)), pnorm(a), tolerance = 2e-12)
    expect_equal(log_bvn_value(a, b, r), log_bvn_value(b, a, r), tolerance = 2e-12)
  }
})

test_that("rare joint events remain finite on the log scale without floors", {
  for (setting in list(c(-8, -7, -.3), c(-35, -30, -.5), c(-12, 3, .8))) {
    a = setting[1]; b = setting[2]; r = setting[3]
    value = log_bvn_value(a, b, r)
    expect_true(is.finite(value))
    expect_lte(value, min(pnorm(a, log.p = TRUE), pnorm(b, log.p = TRUE)))
    expect_equal(value, log_bvn_value(b, a, r), tolerance = 1e-10)
    x = jet(a, 3); y = jet(b, 3); correlation = jet(r, 3)
    x$gradient[1, 1] = y$gradient[1, 2] = correlation$gradient[1, 3] = 1
    gradient = jet_log_bvn(x, y, correlation)$gradient
    numerical = core_difference(function(v) log_bvn_value(v[1], v[2], v[3]), setting,
                                 step = 1e-6)
    expect_lt(max(abs(drop(gradient) - numerical)/(1 + abs(numerical))), 3e-6)
  }
  expect_lt(log_bvn_value(-35, -30, -.5), log(.Machine$double.xmin))
  extreme = log_bvn_value(-35, -30, -.9999)
  expect_true(is.finite(extreme))
  expect_lt(extreme, -1e7)
  expect_equal(extreme, log_bvn_value(-30, -35, -.9999))
})

test_that("tail mode bracket failures use the existing invalid-trial condition", {
  # Captured from a representable, extreme optimizer trial in bootstrap b059.
  arguments = c(a = 7562158621.9907141, b = 56419477043.966141,
                r = -0.99999999997757327)
  fail = function() do.call(log_bvn_tail, as.list(arguments))
  expect_error(fail(), "Bivariate normal tail mode: f\\(\\) values at end points",
               class = "panelselect_integration_error")

  fixture = core_fixture("selection_probit")
  info = parameter_layout(fixture$prepared)
  info$fixed = FALSE
  theta = fixture$theta
  local_mocked_bindings(evaluate_model = function(...) fail(), .package = "PanelSelect")
  evaluator = fit_evaluator(theta, info, fixture$prepared, core_control("agh"), numeric())
  eta = parameter_to_optimizer(theta, info)
  value = evaluator$evaluate(eta, trial = TRUE)
  expect_false(value$valid)
  expect_s3_class(value$error, "panelselect_integration_error")
  expect_true(all(value$loglik == -Inf))
  expect_true(all(value$optimizer_score == 0))
  expect_identical(evaluator$rejected(), 1L)
  # The cached rejected trial remains an explicit error for a final evaluation.
  expect_error(evaluator$evaluate(eta, trial = FALSE), "Bivariate normal tail mode:",
               class = "panelselect_integration_error")
  expect_identical(evaluator$rejected(), 1L)
  eta[1] = eta[1]+0.01
  expect_error(evaluator$evaluate(eta, trial = FALSE), "Bivariate normal tail mode:",
               class = "panelselect_integration_error")
  expect_identical(evaluator$rejected(), 2L)
})
