# Numerical values are checked by the independent fixtures and derivative tests.
# These probes ensure each model file's own functions are the executable path:
# its panel kernel at person-effect nodes and, for the fixed-rule shortcuts,
# its outcome rows.
owned_models = function() {
  original_resolve = resolve_model
  function(...) {
    model = original_resolve(...)
    expect_true(is.function(model$panel), info = model$id)
    model$panel = function(...) stop('panel kernel of ', model$id, ' reached')
    if (!is.null(model$outcome_rows))
      model$outcome_rows = function(...) stop('outcome rows of ', model$id, ' reached')
    model
  }
}

test_that('every model evaluates its own panel kernel at person-effect nodes', {
  local_mocked_bindings(resolve_model = owned_models())
  cases = lapply(core_models, core_fixture)
  for (family in c('probit', 'gaussian', 'pln')) {
    fixture = core_fixture(paste0('selection_', family))
    fixture$prepared$structure = 'cross_sectional'
    fixture$prepared$model_id = c(probit = 'probit_probit', gaussian = 'probit_linear',
      pln = 'probit_pln')[[family]]
    fixture$theta = fixture$theta[!names(fixture$theta) %in%
      c('sd_re_selection', 'sd_re_outcome', 'cor_re')]
    cases[[length(cases)+1L]] = fixture
  }
  for (fixture in cases) {
    control = core_control(nodes = 3L, inner = 5L)
    q = model_dimension(fixture$theta, fixture$prepared, control, numeric())
    z = matrix(.1, fixture$prepared$n_panels, q)
    id = fixture$prepared$model_id
    for (order in 0:3) {
      kernel = make_panel_kernel(fixture$theta, fixture$prepared, control, numeric(),
        order = order, parameters = order > 0L, latent = TRUE)
      # Standalone models evaluate first-order nodes with the shared shortcut,
      # which calls their outcome rows.
      expected = if (!fixture$prepared$has_selection && order <= 1L)
        paste('outcome rows of', id, 'reached') else paste('panel kernel of', id, 'reached')
      expect_error(kernel$evaluate(z), expected, fixed = TRUE)
    }
  }
})

test_that('fixed-rule shortcuts evaluate the model-owned outcome rows', {
  local_mocked_bindings(resolve_model = owned_models())
  control = core_control('gh', nodes = 3L, inner = 5L)
  # Standalone fixed GH averages the outcome rows' node partials before
  # expanding the design.
  for (name in c('probit_re', 'poisson_re', 'pln_re')) {
    fixture = core_fixture(name)
    for (score in c(FALSE, TRUE))
      expect_error(evaluate_model(fixture$theta, fixture$prepared, control, score = score),
        paste('outcome rows of', name, 'reached'), fixed = TRUE)
  }
  # With tau fixed at 0, fixed GH integrates the outcome rows given each
  # selection node.
  for (name in c('selection_probit', 'selection_poisson', 'selection_pln')) {
    fixture = core_fixture(name)
    fixed = if ('cor_error' %in% names(fixture$theta)) c(cor_error = 0) else numeric()
    theta = fixture$theta; theta[names(fixed)] = fixed
    expect_error(evaluate_model(theta, fixture$prepared, control, fixed),
      paste('outcome rows of', fixture$prepared$model_id, 'reached'), fixed = TRUE)
  }
})
