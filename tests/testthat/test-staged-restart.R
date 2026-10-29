# Stage ownership with a real evaluator. Scripted likelihood values let each
# test define the expected argmax or count itself.
staged_restart_fixture = function() {
  theta = c('outcome.(Intercept)' = 0)
  list(theta = theta,
    info = data.frame(name = names(theta), type = 'coefficient', fixed = FALSE),
    prepared = list(n_panels = 1L, n_rows = 2L, family = 'probit',
      model_id = 'probit_re', structure = 'panel', has_selection = FALSE,
      x = matrix(1, 2L, 1L, dimnames = list(NULL, '(Intercept)')),
      w = matrix(numeric(), 2L, 0L), selected = c(TRUE, TRUE)),
    fixed = numeric(), control = fit_control(panel_control(quadrature = 'gl')))
}

test_that('a staged restart begins at the best valid point of the failed stage', {
  fixture = staged_restart_fixture()
  # The best valid trial (-5 at theta 2) is not the last valid trial (-7 at 3).
  values = c('0' = -20, '1' = -10, '2' = -5, '3' = -7)
  starts = list(); stage = 0L
  local_mocked_bindings(
    uses_staged_integration = function(...) TRUE,
    prepare_model_data = function(prepared) prepared,
    prepare_staged_integration = function(theta, prepared, control, fixed,
        probes = NULL, previous = NULL, epoch = 1L) {
      stage <<- epoch; starts[[epoch]] <<- theta
      list(epoch = epoch)
    },
    evaluate_model = function(theta, prepared, control, fixed, score = TRUE) {
      key = format(unname(theta))
      if (!key %in% names(values))
        stop(staged_refinement_condition('Deliberate trial refinement',
          probes = list(list(theta = theta))))
      score = matrix(0, 1L, 1L, dimnames = list(NULL, names(theta)))
      list(loglik = values[[key]], score = score, reference_score = score,
        diagnostics = list(method = 'scripted'))
    },
    optimize_model = function(eta, evaluator, control, options) {
      if (stage == 1L) for (value in 1:4) evaluator$evaluate(stats::setNames(value, names(eta)))
      list(par = eta, counts = c('function' = 1L, gradient = 1L), convergence = 0L, message = NULL)
    },
    validate_staged_integration = function(...) invisible(TRUE),
    .package = 'PanelSelect')
  fit = do.call(fit_staged_model, fixture)
  records = fit$final$diagnostics$staging$records
  expect_length(records, 2L)
  expect_identical(unname(records[[2L]]$start), 2)
  expect_identical(unname(starts[[2L]]), 2)
})

test_that('a failure raised by mesh preparation is terminal, not retried', {
  fixture = staged_restart_fixture()
  calls = 0L
  local_mocked_bindings(prepare_staged_integration = function(...) {
    calls <<- calls + 1L
    stop(staged_refinement_condition('Deliberate tail unresolved at maximum bound'))
  }, .package = 'PanelSelect')
  failure = tryCatch(do.call(fit_staged_model, fixture), error = identity)
  expect_s3_class(failure, 'panelselect_mesh_required')
  expect_identical(calls, 1L)
  expect_identical(failure$diagnostics$stages_used, 1L)
  expect_identical(failure$diagnostics$phase, 'prepare')
  expect_false(grepl('after three total stages', conditionMessage(failure)))
})

test_that('a real split preparation failure at its fixed bound stops after one preparation', {
  # A truncation bound of 2.5 cannot meet the default tail tolerance for any
  # data, so the failure is deterministic. No random numbers are used.
  id = rep(1:10, each = 3L)
  x = cos(1.3*seq_along(id))
  y = floor(exp(.2 + .4*x + 1.5*stats::qnorm(stats::ppoints(10))[id]))
  data = data.frame(id = id, x = x, y = y)
  original = prepare_staged_integration
  calls = 0L
  local_mocked_bindings(prepare_staged_integration = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = 'PanelSelect')
  failure = tryCatch(panel_select(y ~ x, NULL, data, 'id', 'poisson',
    control = panel_control(quadrature = 'gl', quadrature_options = list(bound = 2.5))), error = identity)
  expect_s3_class(failure, 'panelselect_mesh_required')
  expect_identical(calls, 1L)
  expect_identical(failure$diagnostics$stages_used, 1L)
  expect_false(grepl('after three total stages', conditionMessage(failure)))
})

test_that('a successful staged fit does not repeat its endpoint evaluation', {
  # Endpoint accuracy is enforced by the final evaluation itself; a separate
  # validation pass repeats the same deterministic integral.
  id = rep(1:6, each = 3L)
  x = sin(seq_along(id))
  data = data.frame(id = id, x = x, y = c(0, 1, 2, 1, 0, 3, 2, 2, 1, 0, 0, 1, 4, 2, 3, 1, 1, 0))
  calls = 0L
  original = validate_staged_integration
  local_mocked_bindings(validate_staged_integration = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = 'PanelSelect')
  budget = list(seconds_per_panel = 1e4, seconds_total = 1e4)
  for (vcov in c('opg', 'observed')) {
    calls = 0L
    fit = panel_select(y ~ x, NULL, data, 'id', 'poisson',
      control = panel_control(quadrature = 'gl', quadrature_options = budget,
        stationarity_warning = NULL), vcov = vcov)
    expect_identical(calls, 0L)
    phases = vapply(fit$integration$staging$records, function(record) record$phase, character(1L))
    expect_identical(tail(phases, 1L), if (vcov == 'opg') 'endpoint_validation' else 'information')
  }
})
