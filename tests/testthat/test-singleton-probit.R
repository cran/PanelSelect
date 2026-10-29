test_that("singleton probit uses actual fixed coefficients and offsets before initialization", {
  data = data.frame(id = 1:6, x = c(-1, -.4, 0, .3, .7, 1.2), y = rep(0:1, 3), off = 0)
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  expect_error(panel_select(y ~ x, NULL, data, "id", "probit"), "known index is in the span")
  expect_error(panel_select(y ~ x, NULL, data, "id", "probit",
    fixed = c("outcome.(Intercept)" = 0)), "known index is in the span")
  data$off = .4 * data$x
  expect_error(panel_select(y ~ 0 + x + offset(off), NULL, data, "id", "probit"),
    "known index is in the span")
  data$off = -.3
  expect_error(panel_select(y ~ x + offset(off), NULL, data, "id", "probit",
    fixed = c("outcome.(Intercept)" = .3)), "known index is in the span")
  expect_error(panel_select(y ~ 1, NULL, data, "id", "probit",
    fixed = c("outcome.(Intercept)" = 0)), "known index is in the span")
  expect_error(panel_select(y ~ 0, NULL, data, "id", "probit"), "known index is in the span")
})

test_that("independent known indices, fixed SDs and repeated panels remain allowed", {
  data = data.frame(id = 1:6, x = c(-1, -.4, 0, .3, .7, 1.2), y = rep(0:1, 3), off = .4)
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  for (value in c(-.3, .3))
    expect_error(panel_select(y ~ x, NULL, data, "id", "probit",
      fixed = setNames(value, "outcome.(Intercept)")), "Initialization reached")
  expect_error(panel_select(y ~ 0 + x + offset(off), NULL, data, "id", "probit"), "Initialization reached")
  expect_error(panel_select(y ~ 1, NULL, data, "id", "probit",
    fixed = c("outcome.(Intercept)" = .3)), "Initialization reached")
  expect_error(panel_select(y ~ 0 + offset(off), NULL, data, "id", "probit"), "Initialization reached")
  for (value in c(0, .7))
    expect_error(panel_select(y ~ x, NULL, data, "id", "probit", fixed = c(sd_re_outcome = value)),
      "Initialization reached")
  data$id[2] = data$id[1]
  expect_error(panel_select(y ~ x, NULL, data, "id", "probit"), "Initialization reached")
})

test_that("rank and response checks run first and singleton detection uses retained rows", {
  data = data.frame(id = rep(1:3, each = 2), x = rep(-1:1, each = 2),
                     y = c(0, NA, 1, NA, 0, NA))
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  expect_error(panel_select(y ~ x, NULL, data, "id", "probit", na_action = "omit"), "known index is in the span")
  data$id = seq_len(nrow(data)); data$y = rep(0:1, 3); data$duplicate = data$x
  expect_error(panel_select(y ~ x + duplicate, NULL, data, "id", "probit"),
    "free outcome coefficients are not identified")
  data$y = 1
  expect_error(panel_select(y ~ x, NULL, data, "id", "probit"), "no finite maximum-likelihood")
})

test_that("the approved backward-error rule is scale invariant near the span", {
  q1 = rep(1, 8) / sqrt(8)
  q2 = rep(c(1, -1), 4) / sqrt(8)
  q3 = rep(c(1, 1, -1, -1), 2) / sqrt(8)
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  for (scale in list(c(1, 1, 1), c(1e-80, 1e80, 1e-100), c(1e80, 1e-80, 1e100))) {
    for (delta in c(0, 1e-13, 1e-12, 1e-8)) {
      data = data.frame(id = 1:8, y = rep(0:1, 4), a = q1 * scale[1], b = q2 * scale[2],
                         off = (.6 * q1 + .8 * q2 + delta * q3) * scale[3])
      expected = if (delta <= 1e-13) "known index is in the span" else "Initialization reached"
      expect_error(panel_select(y ~ 0 + a + b + offset(off), NULL, data, "id", "probit"), expected)
    }
  }
  # Cancellation in a nearly dependent design requires the coefficient-norm
  # denominator; the unscaled projection residual alone would miss this case.
  data = data.frame(id = 1:8, y = rep(0:1, 4), a = q1, b = q1 + 1e-5 * q2)
  data$off = data$b - data$a
  expect_error(panel_select(y ~ 0 + a + b + offset(off), NULL, data, "id", "probit"),
    "known index is in the span")
})

test_that("direct normalization handles extreme finite known indices", {
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  data = data.frame(id = 1:8, y = rep(0:1, 4), x = rep(c(1, -1), 4), off = 0)
  for (scale in c(1e308, 1e-308)) {
    data$off = rep(scale, 8)
    expect_error(panel_select(y ~ 0 + x + offset(off), NULL, data, "id", "probit"), "Initialization reached")
    data$off = data$x * scale
    expect_error(panel_select(y ~ 0 + x + offset(off), NULL, data, "id", "probit"),
      "known index is in the span")
  }
  data$x = 2; data$off = 0
  expect_error(panel_select(y ~ 0 + x, NULL, data, "id", "probit", fixed = c(outcome.x = 1e308)),
      "fixed probit index is not finite")
})

test_that("large cancelling representations do not overflow the backward residual", {
  # This QR-full-rank design has representation coefficients around 1e190.
  # Squaring those coefficients would return NA instead of a policy decision.
  x = diag(400)
  x[upper.tri(x)] = -2
  x = rbind(x, 0)
  known = c(rep(0, 399), 1, .2)
  expect_equal(qr(x)$rank, ncol(x))
  expect_true(in_design_span(x, known))
})
