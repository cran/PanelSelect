test_that('cached integration respects changed controls and leaves inputs unchanged', {
  fixture = core_fixture('selection_pln')
  control = core_control('gh', nodes = 3L, inner = 5L)
  before = serialize(list(fixture, control), NULL)
  prepared = prepare_integration(fixture$prepared, control)
  expect_identical(serialize(list(fixture, control), NULL), before)
  expect_identical(prepare_integration(prepared, control), prepared)
  for (method in c('gh', 'agh')) for (reduction in c(TRUE, FALSE)) {
    changed = control
    changed$quadrature = method
    changed$exact_reductions = reduction
    changed$mode_tol = 2e-10
    cached = evaluate_model(fixture$theta, prepared, changed)
    fresh = evaluate_model(fixture$theta, fixture$prepared, changed)
    expect_equal(cached, fresh, tolerance = 1e-12)
  }
})

test_that('node blocks preserve likelihoods, scores and posterior ordering', {
  one_node = function(...) 1L
  for (model in core_models) for (method in c('gh', 'agh')) {
    fixture = core_fixture(model)
    control = core_control(method, nodes = 3L, inner = 5L, exact_reductions = FALSE)
    blocked = evaluate_model(fixture$theta, fixture$prepared, control)
    serial = with_mocked_bindings(
      evaluate_model(fixture$theta, fixture$prepared, control),
      quadrature_block_size = one_node)
    expect_equal(blocked$loglik, serial$loglik, tolerance = 2e-12, info = model)
    expect_equal(blocked$score, serial$score, tolerance = 2e-11, info = model)
    if (!model %in% c('selection_gaussian', 'selection_pln')) next
    blocked = posterior_rule(fixture$theta, fixture$prepared, control)
    serial = with_mocked_bindings(
      posterior_rule(fixture$theta, fixture$prepared, control),
      quadrature_block_size = one_node)
    expect_equal(blocked, serial, tolerance = 2e-11, info = model)
  }
})

test_that('component memo keys follow equation dependencies and evaluation settings', {
  fixture = core_fixture('selection_pln')
  control = core_control('gh', nodes = 3L, inner = 5L)
  prepared = prepare_integration(fixture$prepared, control)
  selection = prepared$integration_cache$selection
  outcome = prepared$integration_cache$outcome
  calls = 0L
  local_mocked_bindings(integrate_prepared = function(theta, prepared, control, fixed, score) {
    calls <<- calls + 1L
    list(loglik = rep(sum(theta), prepared$n_panels),
      score = matrix(if (score) rep(theta, each = prepared$n_panels) else 0,
        prepared$n_panels, length(theta)))
  })
  evaluate = function(theta = fixture$theta, component = selection, fixed = numeric(),
                      score = TRUE, settings = control)
    integrate_component(theta, component, settings, fixed, score)

  first = evaluate()
  changed = fixture$theta; changed['outcome.x'] = changed['outcome.x'] + .1
  expect_identical(evaluate(changed), first)
  expect_identical(calls, 1L)
  expect_identical(colnames(first$score), names(fixture$theta))
  irrelevant = !startsWith(names(fixture$theta), 'selection.') &
    names(fixture$theta) != 'sd_re_selection'
  expect_true(all(first$score[, irrelevant] == 0))
  changed['selection.w'] = changed['selection.w'] + .1
  expect_false(identical(evaluate(changed), first))
  expect_identical(calls, 2L)

  initial_outcome = evaluate(component = outcome)
  changed = fixture$theta; changed['selection.w'] = changed['selection.w'] + .1
  expect_identical(evaluate(changed, outcome), initial_outcome)
  expect_identical(calls, 3L)
  changed['sd_error'] = changed['sd_error'] + .1
  expect_false(identical(evaluate(changed, outcome), initial_outcome))
  expect_identical(calls, 4L)

  expect_true(all(evaluate(score = FALSE)$score == 0))
  expect_identical(evaluate(), first)
  expect_identical(calls, 6L)
  fixed = fixture$theta['sd_re_selection']
  expect_identical(evaluate(fixed = fixed), first)
  expect_identical(evaluate(fixed = fixed), first)
  expect_identical(calls, 7L)

  changed_control = control; changed_control$nodes = 4L
  fresh = prepare_integration(prepared, changed_control)$integration_cache$selection
  expect_false(identical(fresh$integration_cache$memo, selection$integration_cache$memo))
  expect_identical(evaluate(component = fresh, fixed = fixed, settings = changed_control), first)
  expect_identical(calls, 8L)
  other = prepare_integration(fixture$prepared, control)$integration_cache$selection
  expect_identical(evaluate(component = other, fixed = fixed), first)
  expect_identical(evaluate(fixed = fixed), first)
  expect_identical(calls, 9L)
  changed = fixture$theta; changed['sd_re_selection'] = changed['sd_re_selection'] * .9
  expect_false(identical(evaluate(changed), first))
  expect_identical(evaluate(component = other, fixed = fixed), first)
  expect_identical(calls, 10L)
})

test_that('failed component trials do not replace a cached successful value', {
  fixture = core_fixture('selection_probit')
  control = core_control('gh', nodes = 3L)
  prepared = prepare_integration(fixture$prepared, control)$integration_cache$selection
  calls = 0L
  local_mocked_bindings(integrate_prepared = function(theta, prepared, control, fixed, score) {
    calls <<- calls + 1L
    if (theta[['selection.w']] > 100) integration_error('Invalid component trial')
    list(loglik = 1, score = matrix(1, prepared$n_panels, length(theta)))
  })
  evaluate = function(theta) integrate_component(theta, prepared, control, numeric(), TRUE)
  valid = evaluate(fixture$theta)
  key = prepared$integration_cache$memo$key
  value = prepared$integration_cache$memo$value
  invalid = fixture$theta; invalid['selection.w'] = 101
  expect_error(evaluate(invalid), 'Invalid component trial',
               class = 'panelselect_integration_error')
  expect_identical(prepared$integration_cache$memo$key, key)
  expect_identical(prepared$integration_cache$memo$value, value)
  expect_identical(evaluate(fixture$theta), valid)
  expect_identical(calls, 2L)
  expect_error(evaluate(invalid), 'Invalid component trial',
               class = 'panelselect_integration_error')
  expect_identical(calls, 3L)
})

test_that('fixed independent panels share selection and preserve outcome partitions', {
  for (model in c('selection_gaussian', 'selection_pln')) {
    fixture = core_fixture(model); prepared = fixture$prepared
    prepared$group = c(1L, 2L, 2L); prepared$group_ptr = c(0L, 1L, 3L)
    prepared$id_values = c('unselected', 'observed'); prepared$n_panels = 2L
    prepared$d = c(0, 1, 1); prepared$selected = prepared$d == 1
    prepared$y = c(NA_real_, fixture$prepared$y[c(1L, 3L)])
    fixed = c(cor_re = 0, cor_error = 0)
    theta = fixture$theta; theta[names(fixed)] = fixed
    control = core_control('gh', nodes = 3L, inner = 5L)
    selection = integrate_prepared(theta, selection_component(prepared), control, fixed)
    part = subset_panels(prepared, 2L)
    if (model == 'selection_gaussian') {
      outcome = linear_outcome_marginal(theta, part, conditional = FALSE)
      outcome = list(loglik = outcome$value, score = outcome$gradient)
    } else {
      part$has_selection = FALSE; part$model_id = 'pln_re'; part$active = part$d == 1
      outcome = integrate_prepared(theta, part, control, fixed)
    }
    expected = selection
    expected$loglik[2L] = expected$loglik[2L] + outcome$loglik
    expected$score[2L, ] = expected$score[2L, ] + outcome$score
    answer = evaluate_model(theta, prepared, control, fixed)
    expect_equal(answer$loglik, expected$loglik, tolerance = 2e-12, info = model)
    expect_equal(answer$score, expected$score, tolerance = 2e-11, info = model)
    irrelevant = !startsWith(names(theta), 'selection.') & names(theta) != 'sd_re_selection'
    expect_true(all(answer$score[1L, irrelevant] == 0))
    expect_identical(answer$diagnostics$selection_only_panels, 'unselected')
    expect_identical(answer$diagnostics$observed_panels, 'observed')
    expect_identical(answer$diagnostics$selection_only, selection$diagnostics)
    if (model == 'selection_gaussian')
      expect_identical(answer$diagnostics$observed, selection$diagnostics)
    else {
      expect_true(answer$diagnostics$observed$factorized)
      expect_identical(answer$diagnostics$observed$selection, selection$diagnostics)
      expect_identical(answer$diagnostics$observed$outcome, outcome$diagnostics)
    }
    prepared$d[] = 0; prepared$selected[] = FALSE; prepared$y[] = NA_real_
    selection = integrate_prepared(theta, selection_component(prepared), control, fixed)
    answer = evaluate_model(theta, prepared, control, fixed)
    expect_equal(answer$loglik, selection$loglik, tolerance = 2e-12)
    expect_equal(answer$score, selection$score, tolerance = 2e-11)
    expect_identical(answer$diagnostics$selection_only_panels, prepared$id_values)
    expect_null(answer$diagnostics$observed)
    expect_identical(answer$diagnostics$observed_panels, character())
  }
})

test_that('Gaussian outcome memo follows outcome coordinates and gradient layout', {
  fixture = core_fixture('selection_gaussian')
  prepared = prepare_integration(fixture$prepared, core_control())
  original_sum = group_sum; calls = 0L
  local_mocked_bindings(group_sum = function(...) {
    calls <<- calls + 1L
    original_sum(...)
  })
  theta = fixture$theta
  marginal = linear_outcome_marginal(theta, prepared, conditional = FALSE)
  before = calls
  changed = theta
  irrelevant = startsWith(names(theta), 'selection.') |
    names(theta) %in% c('sd_re_selection', 'cor_re', 'cor_error')
  changed[irrelevant] = changed[irrelevant] + .01
  expect_identical(linear_outcome_marginal(changed, prepared, conditional = FALSE), marginal)
  expect_identical(calls, before)
  for (name in c('outcome.x', 'sd_re_outcome', 'sd_error')) {
    changed = theta; changed[name] = changed[name] + .01
    before = calls
    next_value = linear_outcome_marginal(changed, prepared, conditional = FALSE)
    expect_gt(calls, before)
    expect_false(identical(next_value$value, marginal$value))
  }
  linear_outcome_marginal(theta, prepared, conditional = FALSE)
  before = calls
  conditional = linear_outcome_marginal(theta, prepared, conditional = TRUE)
  expect_gt(calls, before)
  expect_equal(conditional[c('value', 'gradient')], marginal)
  expect_true(all(c('mean', 'variance', 'scale', 'mean_gradient',
    'variance_gradient', 'scale_gradient') %in% names(conditional)))
  before = calls
  reordered = linear_outcome_marginal(rev(theta), prepared, conditional = TRUE)
  expect_gt(calls, before)
  for (field in c('gradient', 'mean_gradient', 'variance_gradient', 'scale_gradient'))
    expect_equal(reordered[[field]], conditional[[field]][, rev(seq_along(theta)), drop = FALSE])
  expect_equal(reordered$value, conditional$value)
})
