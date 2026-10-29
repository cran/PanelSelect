test_that("constant binary responses and zero counts expose a monotone likelihood path", {
  for (model in setdiff(core_models, "selection_gaussian")) {
    fixture = core_fixture(model)
    prepared = fixture$prepared
    for (value in if (prepared$family == "probit") 0:1 else 0) {
      prepared$y[prepared$selected] = value
      info = parameter_layout(prepared)
      info$fixed = info$name != "outcome.(Intercept)"
      expect_error(validate_design_rank(prepared, info), "outcome response has no finite")
      likelihood = vapply(0:2, function(step) {
        theta = fixture$theta
        theta["outcome.(Intercept)"] = (2*value-1)*step
        sum(evaluate_model(theta, prepared, core_control(nodes = 3, inner = 5),
                           score = FALSE)$loglik)
      }, numeric(1))
      expect_true(all(diff(likelihood) > 0))
    }
  }
})

test_that("constant participation exposes the same path for every selected model", {
  for (model in grep("^selection_", core_models, value = TRUE)) {
    fixture = core_fixture(model)
    for (value in 0:1) {
      prepared = fixture$prepared
      prepared$d[] = value
      prepared$selected[] = value == 1
      prepared$y[] = 0
      info = parameter_layout(prepared)
      info$fixed = info$name != "selection.(Intercept)"
      expect_error(validate_design_rank(prepared, info), "selection response has no finite")
      likelihood = vapply(0:2, function(step) {
        theta = fixture$theta
        theta["selection.(Intercept)"] = (2*value-1)*step
        sum(evaluate_model(theta, prepared, core_control(nodes = 3, inner = 5),
                           score = FALSE)$loglik)
      }, numeric(1))
      expect_true(all(diff(likelihood) > 0))
    }
  }
})

test_that("the guard detects a constant span without requiring an intercept column", {
  data = data.frame(id = rep(1:2, each = 2), y = 0, group = factor(c(1, 2, 1, 2)),
                     known = c(-.2, .1, .4, -.3))
  expect_error(panel_select(y ~ 0 + group + offset(known), NULL, data, "id", "poisson",
    fixed = c(sd_re_outcome = 0)), "free constant predictor direction")
  data$y[] = 2
  prepared = prepare_data(y ~ 0 + group, NULL, data, "id", "poisson")
  info = parameter_layout(prepared)
  expect_silent(validate_design_rank(prepared, info))
})

test_that("fixed coefficients and known offsets preserve constant-response likelihoods", {
  for (model in core_models) {
    fixture = core_fixture(model)
    data = data.frame(id = 1, y = 0, d = 0, x = c(-.4, .3, .7), w = c(.5, -.3, .1))
    selected = fixture$prepared$has_selection
    for (value in if (selected) 0:1 else 0) {
      data$d[] = value
      fit = panel_select(y ~ x,
        if (selected) d ~ w else NULL, data, "id", fixture$prepared$family, fixed = fixture$theta,
        control = panel_control(quadrature = "gh", nodes = 3, observation_nodes = 5))
      expect_true(is.finite(as.numeric(logLik(fit))))
      expect_equal(attr(logLik(fit), "df"), 0L)
    }
  }
  # One panel makes the OPG-scaled score statistic identically one, whatever the
  # gradient. With every scale fixed at zero, separate panels leave the
  # likelihood unchanged and keep that diagnostic informative.
  data = data.frame(id = 1:3, y = 0, d = 1, x = c(-1, 0, 1), known = c(.2, -.2, .1))
  for (family in c("probit", "poisson", "pln")) {
    fixed = c("outcome.(Intercept)" = 0, sd_re_outcome = 0,
               if (family == "pln") c(sd_error = 0))
    fit = panel_select(y ~ x + offset(known), NULL, data, "id", family, fixed = fixed)
    expect_true(fit$optimizer$converged)
    expect_true(is.finite(coef(fit)[["outcome.x"]]))
  }
  fit = panel_select(y ~ 1,
    d ~ x + offset(known), data, "id", "gaussian",
    fixed = c("outcome.(Intercept)" = 0, "selection.(Intercept)" = 0,
              sd_re_outcome = 0, sd_re_selection = 0, sd_error = 1, cor_error = 0))
  expect_true(fit$optimizer$converged)
  expect_true(is.finite(coef(fit)[["selection.x"]]))
})
