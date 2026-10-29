# Coverage for stored observed information, prediction intervals, quadrature
# diagnostics with an independent GL rule, terms()/na.action() and control printing.

capabilities_data = function(seed = 3201, panels = 18L, waves = 3L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    z = stats::rnorm(n))
  u = stats::rnorm(panels)
  data$d = as.integer(.3 + .5*data$x + .7*data$z + .6*u[data$id] + stats::rnorm(n) > 0)
  data$count = stats::rpois(n, exp(.2 + .4*data$x + .5*u[data$id]))
  data$binary = as.integer(.1 + .5*data$x + .5*u[data$id] + stats::rnorm(n) > 0)
  data$linear = .2 + .5*data$x + .5*u[data$id] + stats::rnorm(n)
  for (name in c("count", "binary", "linear")) data[[name]][data$d == 0] = NA
  data
}

capabilities_fit = function(data, family = "poisson", vcov = "opg", ...) {
  outcome = stats::as.formula(paste(switch(family, poisson = "count", probit = "binary",
    gaussian = "linear"), "~ x"))
  panel_select(outcome, d ~ x + z, data, "id", family, fixed = c(cor_re = .3, sd_re_outcome = .5),
    control = panel_control(nodes = 4L, reltol = 1e-10, ...), init = "pooled", vcov = vcov)
}

test_that("add_information() stores what an observed fit stores, once", {
  data = capabilities_data()
  fit = capabilities_fit(data)
  expect_false(fit$optimizer$information_computed)
  informed = add_information(fit)
  expect_true(informed$optimizer$information_computed)
  expect_identical(vcov(informed, "observed"), vcov(fit, "observed"))
  expect_identical(vcov(informed, "sandwich"), vcov(fit, "sandwich"))
  expect_identical(vcov(informed), vcov(fit))
  expect_identical(informed$control$vcov, "opg")
  expect_identical(add_information(informed), informed)
  expect_identical(check_stationarity(fit)$status, "unavailable")
  expect_identical(check_stationarity(informed)$status, "available")
  expect_identical(predict(informed, se.fit = TRUE, vcov = "observed"),
    predict(fit, se.fit = TRUE, vcov = "observed"))
  restored = unserialize(serialize(informed, NULL))
  expect_identical(vcov(restored, "observed"), vcov(informed, "observed"))
  # A fit computed with vcov = "observed" already stores the same information.
  observed = capabilities_fit(data, vcov = "observed")
  expect_identical(add_information(observed), observed)
  expect_equal(informed$optimizer$information, observed$optimizer$information, tolerance = 1e-12)
  # An older object without the newer fields still gains its information.
  old = fit
  old$optimizer$evaluations = NULL
  old$prepared$design$outcome$assign = old$prepared$design$selection$assign = NULL
  expect_identical(vcov(add_information(old), "observed"), vcov(informed, "observed"))
})

test_that("a staged fit whose mesh cannot serve the stencil records why", {
  data = capabilities_data(panels = 6L)
  fit = panel_select(count ~ x, d ~ x, data, "id", "poisson", fixed = c(cor_re = .3, sd_re_selection = .5,
    sd_re_outcome = .4, "selection.(Intercept)" = .1, selection.x = .2, "outcome.(Intercept)" = .1,
    outcome.x = .3), control = panel_control(quadrature = "gl"))
  fit$optimizer$information_computed = FALSE
  local_mocked_bindings(validate_staged_integration = function(...)
    stop(staged_refinement_condition("stencil needs a finer mesh")), .package = "PanelSelect")
  result = with_warnings(add_information(fit))
  expect_match(result$warnings, "requires a refined integration mesh", all = TRUE)
  expect_identical(result$value$optimizer$information_error, "stencil needs a finer mesh")
  expect_true(result$value$optimizer$information_computed)
})

test_that("confidence intervals use the probit, log and identity links", {
  data = capabilities_data()
  z = stats::qnorm(.975)
  targets = data[1:6, ]
  count = capabilities_fit(data, "poisson")
  value = predict(count, targets, interval = "confidence", se.fit = TRUE)
  m = value$fit[, "fit"]; se = value$se.fit
  expect_identical(colnames(value$fit), c("fit", "lwr", "upr"))
  expect_equal(unname(value$fit[, "lwr"]), unname(exp(log(m) - z*se/m)), tolerance = 1e-12)
  expect_equal(unname(value$fit[, "upr"]), unname(exp(log(m) + z*se/m)), tolerance = 1e-12)
  expect_identical(m, predict(count, targets))
  p = predict(count, targets, type = "selection_probability", interval = "confidence", level = .9, se.fit = TRUE)
  q = stats::qnorm(p$fit[, "fit"]); half = stats::qnorm(.95)*p$se.fit/stats::dnorm(q)
  expect_equal(unname(p$fit[, "lwr"]), unname(stats::pnorm(q - half)), tolerance = 1e-12)
  expect_equal(unname(p$fit[, "upr"]), unname(stats::pnorm(q + half)), tolerance = 1e-12)
  linear = capabilities_fit(data, "gaussian")
  value = predict(linear, targets, interval = "confidence", se.fit = TRUE)
  expect_equal(unname(value$fit[, "upr"] - value$fit[, "fit"]), unname(z*value$se.fit), tolerance = 1e-12)
  link = predict(count, targets, type = "link", interval = "confidence", se.fit = TRUE)
  expect_equal(unname(link$fit[, "lwr"]), unname(link$fit[, "fit"] - z*link$se.fit), tolerance = 1e-12)
  # Aggregates and history conditioning return one bounded row per target.
  average = predict(count, targets, aggregate = TRUE, interval = "confidence")
  expect_identical(dim(average), c(1L, 3L))
  expect_equal(unname(average[, "fit"]), predict(count, targets, aggregate = TRUE), tolerance = 1e-15)
  history = predict(count, targets, conditioning = "history", history = data, interval = "confidence")
  expect_true(all(history[, "lwr"] <= history[, "fit"] & history[, "fit"] <= history[, "upr"]))
  expect_error(predict(count, targets, interval = "prediction"), "not available")
  expect_error(predict(count, targets, interval = "confidence", level = 1), "strictly between")
})

test_that("excluded rows are padded in interval matrices", {
  data = capabilities_data()
  data$x[c(2L, 5L)] = NA
  fit = panel_select(count ~ x, d ~ x + z, data, "id", "poisson",
    fixed = c(cor_re = .3, sd_re_outcome = .5, sd_re_selection = .6),
    control = panel_control(nodes = 4L, reltol = 1e-10), na_action = "exclude", init = "pooled")
  value = predict(fit, interval = "confidence")
  expect_identical(nrow(value), nrow(data))
  expect_true(all(is.na(value[c(2L, 5L), ])))
  expect_false(anyNA(value[-c(2L, 5L), ]))
})

test_that("check_quadrature() returns one classed schema on every route", {
  data = capabilities_data()
  columns = c("rule", "nodes", "reference_nodes", "observation_nodes", "observation_quadrature",
    "mesh_changed",
    "loglik_change", "max_panel_change", "max_score_change", "max_panel_score_change",
    "integration_failed", "within_tolerance", "message")
  gh = check_quadrature(capabilities_fit(data))
  expect_s3_class(gh, "panelselect_quadrature_check")
  expect_identical(names(gh$comparison), columns)
  expect_identical(gh$compare, "higher_order")
  expect_output(print(gh),
    "Quadrature check \\(higher_order\\) of: Fixed Gauss-Hermite \\(gh\\), 4 x 4 nodes")
  theta = c("selection.(Intercept)" = .1, selection.x = .2, selection.z = .3,
    "outcome.(Intercept)" = .1, outcome.x = .3, sd_re_selection = .5, sd_re_outcome = .4,
    cor_re = .3, cor_error = 0)
  split = panel_select(binary ~ x, d ~ x + z, data, "id", "probit", fixed = theta,
    control = panel_control(quadrature = "gl"))
  expect_identical(names(check_quadrature(split, nodes = 30L)$comparison), columns)
  small = capabilities_data(panels = 5L)
  staged = panel_select(count ~ x, d ~ x + z, small, "id", "poisson",
    fixed = theta[setdiff(names(theta), "cor_error")], control = panel_control(quadrature = "gl"))
  refined = check_quadrature(staged)
  expect_identical(names(refined$comparison), columns)
  expect_type(refined$comparison$mesh_changed, "logical")
  expect_error(check_quadrature(staged, quadrature_options = list(bound = 9)), "compare = 'gl'")
})

test_that("compare = 'gl' matches the frozen refined references", {
  for (model in c("probit_re", "poisson_re", "pln_re", "selection_probit", "selection_poisson",
                  "selection_gaussian")) {
    fixture = core_fixture(model)
    reference = core_reference(fixture, core_references$refined[[model]])
    panel = fixture$panel
    data = data.frame(id = 1L, x = panel$x[, 2L], w = panel$w[, 2L], y = panel$y, d = panel$d)
    selection = if (fixture$prepared$has_selection) d ~ w
    fit = panel_select(y ~ x, selection, data, "id", fixture$prepared$family,
      fixed = fixture$theta, control = panel_control(quadrature = "gh", nodes = 7L, observation_nodes = 13L))
    check = check_quadrature(fit, compare = "gl")
    expect_identical(check$compare, "gl")
    expect_identical(check$comparison$rule, "gl")
    gl = fit$loglik + check$comparison$loglik_change
    expect_lt(abs(gl - reference$loglik), 1e-7, label = paste(model, "GL log error"))
  }
  data = capabilities_data()
  cross = panel_select(count ~ x, d ~ x + z, data, NULL, "linear", fixed = c(
    "selection.(Intercept)" = .1, selection.x = .2, selection.z = .3, "outcome.(Intercept)" = .1,
    outcome.x = .3, sd_error = 1, cor_error = 0))
  expect_error(check_quadrature(cross, compare = "gl"), "no numerical integral")
})

test_that("terms(), na.action() and the control print expose stored structure", {
  data = capabilities_data()
  data$x[3L] = NA
  fit = panel_select(count ~ x, d ~ x + z, data, "id", "poisson", fixed = c(cor_re = .3, sd_re_outcome = .5),
    control = panel_control(nodes = 4L, reltol = 1e-10), na_action = "omit", init = "pooled")
  expect_identical(terms(fit, "selection"), attr(model.frame(fit, "selection"), "terms"))
  expect_identical(terms(fit), attr(model.frame(fit), "terms"))
  expect_identical(na.action(fit), attr(model.frame(fit), "na.action"))
  expect_identical(unname(c(na.action(fit))), 3L)
  expect_null(na.action(capabilities_fit(capabilities_data())))
  text = capture.output(print(panel_control()))
  expect_identical(text[1:2], c("PanelSelect control",
    "Integration: fixed Gauss-Hermite (gh), 10 nodes"))
  expect_false(any(startsWith(text, "Fit settings")))
  expect_true("Fit settings: na_action omit, init pooled, vcov opg" %in% capture.output(print(fit$control)))
  expect_true("Integration: Gauss-Legendre (gl), order 20 per segment (reference order 10, bound 12)" %in%
    capture.output(print(panel_control(quadrature = "gl"))))
})
