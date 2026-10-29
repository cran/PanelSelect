# Likelihoods, scores and predictions of PanelSelect 1.1.2 (fixtures/panelselect-1.1.2.R)
# are independent references for the native engine.
panelselect112_prepared = function(case) {
  prepare_data(case$outcome, case$selection, case$data, "id", case$family)
}

test_that("binary and Gaussian selection match PanelSelect 1.1.2 finite rules", {
  fixture = fixture_names(dget(test_path("fixtures", "panelselect-1.1.2.R")))
  for (case in Filter(function(case) !case$tail, fixture$cases)) {
    prepared = panelselect112_prepared(case)
    # General Gaussian conditioning changes the finite grid. At fixed tau=0,
    # both versions use the same exact outcome marginal and conditional 1D rule.
    reduction = case$family != "gaussian" || "cor_error" %in% names(case$fixed)
    control = panel_control(quadrature = "gh", nodes = case$H, exact_reductions = reduction)
    value = evaluate_model(case$theta, prepared, control, case$fixed)
    free = !names(case$theta) %in% names(case$fixed)
    expect_identical(as.character(prepared$id_values), case$reference$panel_ids)
    expect_equal(unname(value$loglik), unname(case$reference$loglik),
                 tolerance = 1e-10, info = case$id)
    expect_equal(unname(value$score[, free, drop = FALSE]),
                 unname(case$reference$score[, free, drop = FALSE]),
                 tolerance = 5e-9, info = case$id)
  }
})

test_that("public predictions and uncertainty match PanelSelect 1.1.2 helpers", {
  fixture = fixture_names(dget(test_path("fixtures", "panelselect-1.1.2.R")))
  for (case in Filter(function(case) !case$tail, fixture$cases)) {
    object = structure(list(coefficients = case$theta, data = case$data,
      prepared = panelselect112_prepared(case)), class = "panelselect")
    data = case$data[, setdiff(names(case$data), c("y", "d")), drop = FALSE]
    reference = case$reference$prediction
    # The frozen reference keeps the 1.1.2 target names.
    for (key in c("response", "selection")) {
      type = c(response = "response", selection = "selection_probability")[[key]]
      result = predict(object, data, type = type, se.fit = TRUE, vcov = case$covariance)
      expect_equal(unname(result$fit), reference[[key]], tolerance = 1e-13, info = case$id)
      expect_equal(unname(result$se.fit), reference[[paste0(key, "_se")]],
                   tolerance = 1e-13, info = case$id)
      expect_equal(unname(result$gradient), unname(reference[[paste0(key, "_gradient")]]),
                   tolerance = 1e-13, info = case$id)
    }
    average = predict(object, data, se.fit = TRUE, aggregate = TRUE, vcov = case$covariance)
    expect_equal(c(unname(average$fit), unname(average$se.fit)), reference$average,
                 tolerance = 1e-13, info = case$id)
  }
})

test_that("the legacy probability floor is replaced by the analytic tail likelihood", {
  fixture = fixture_names(dget(test_path("fixtures", "panelselect-1.1.2.R")))
  case = fixture$cases[["probit-legacy-floor"]]
  index = case$theta[c("selection.(Intercept)", "outcome.(Intercept)")]
  scales = case$theta[c("sd_re_selection", "sd_re_outcome")]
  z = index / sqrt(1 + scales^2)
  reference = sum(pnorm(z, log.p = TRUE))
  mills = exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))
  gradient = stats::setNames(numeric(length(case$theta)), names(case$theta))
  gradient[names(index)] = mills / sqrt(1 + scales^2)
  gradient[names(scales)] = -mills * index * scales / (1 + scales^2)^(3/2)
  value = evaluate_model(case$theta, panelselect112_prepared(case),
    panel_control(quadrature = "gh", nodes = case$H), case$fixed)
  expect_equal(unname(case$reference$loglik), log(1e-100), tolerance = 1e-12)
  expect_equal(unname(case$reference$score), matrix(0, 1L, length(case$theta)))
  expect_gt(abs(sum(case$reference$loglik) - reference), 100)
  expect_equal(sum(value$loglik), reference, tolerance = 1e-10)
  expect_equal(unname(colSums(value$score)), unname(gradient), tolerance = 2e-8)
})
