# Likelihoods, derivatives and predictions executed in PanelCount 2.0.1
# (fixtures/panelcount-2.0.1.R) are independent references for the native engine.
panelcount_prepared = function(case) {
  outcome = if (case$selected) observed ~ x + offset(offset_x) else
    if (case$family == "probit") binary ~ x else count ~ x
  selection = if (case$selected) d ~ w + offset(offset_w) else NULL
  prepare_data(outcome, selection, case$data, "id", case$family)
}

test_that("five native models agree with actual PanelCount fixed-rule likelihoods", {
  fixture = fixture_names(dget(test_path("fixtures", "panelcount-2.0.1.R")))
  expect_gt(length(fixture$cases), 0L)
  for (case in fixture$cases) {
    prepared = panelcount_prepared(case)
    control = panel_control(quadrature = "gh", nodes = case$H, observation_nodes = case$H,
                             exact_reductions = FALSE)
    result = evaluate_model(case$theta, prepared, control, fixed = numeric(), score = TRUE)
    score = colSums(result$score)
    if (!case$correlated) {
      expect_lt(abs(sum(result$loglik) - case$reference$loglik), 1e-11, label = case$id)
      # Correlation perturbations move the old/new tensor grids differently even
      # at zero. These are shared non-correlation finite-objective derivatives.
      comparable = !startsWith(names(case$theta), "cor_")
      expect_lt(max(abs(score[comparable] - case$reference$fd_natural[comparable, 2L])), 3e-7,
        label = case$id)
    } else {
      # Different tensor orientations approach the same integral. The source-bound
      # 12/20/30-node ladder justifies these approximation tolerances at order30.
      expect_lt(abs(sum(result$loglik) - case$reference$loglik), 5e-7, label = case$id)
      expect_lt(max(abs(score - case$reference$raw_natural_score)), 3e-6, label = case$id)
    }
  }
})

test_that("public population predictions agree with executed PanelCount helpers", {
  fixture = fixture_names(dget(test_path("fixtures", "panelcount-2.0.1.R")))
  selected = Filter(function(x) x$selected, fixture$cases)
  expect_gt(length(selected), 0L)
  for (case in selected) {
    # Fixed-parameter prediction objects avoid optimizer and covariance effects.
    object = structure(list(coefficients = case$theta,
      prepared = panelcount_prepared(case), data = case$data), class = "panelselect")
    for (pair in list(c("selection_probability", "prob"), c("response", "outcome"), c("joint", "actual_outcome")))
      expect_equal(unname(predict(object, newdata = case$data, type = pair[1L])),
        unname(case$reference$prediction[[pair[2L]]]), tolerance = 1e-13, info = case$id)
  }
})
