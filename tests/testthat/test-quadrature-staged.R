# Staged GL person-effect rule with the observation rule `inner`.
staged_test_control = function(inner = 'gh') {
  staged_test_budget(panel_control(quadrature = 'gl', nodes = 15L,
    quadrature_options = list(reference_nodes = 20L),
    observation_quadrature = inner, observation_nodes = 15L,
    observation_options = if (inner == 'gl') list(reference_nodes = 20L) else list()))
}

# This fixture is local to this test file; test order must not supply data.
pln_inner_fixture = function(model = 'pln_re') {
  selected = model != 'pln_re'; panel = model != 'probit_pln'
  theta = c('outcome.(Intercept)' = .25)
  if (selected) theta = c('selection.(Intercept)' = -.3, theta)
  if (panel) theta = c(theta, sd_re_outcome = .7)
  if (selected && panel) theta = c(theta, sd_re_selection = .6, cor_re = .35)
  theta = c(theta, sd_error = .5)
  if (selected) theta = c(theta, cor_error = -.25)
  one = matrix(1, 1L, 1L, dimnames = list(NULL, '(Intercept)'))
  prepared = list(model_id = model, structure = if (panel) 'panel' else 'cross_sectional',
    family = 'pln', has_selection = selected, x = one,
    w = if (selected) one else matrix(numeric(), 1L, 0L), y = 10, d = 1,
    offset_outcome = 0, offset_selection = 0, group = 1L, group_ptr = c(0L,1L),
    id_values = 'one', n_rows = 1L, n_panels = 1L, selected = TRUE)
  list(theta = theta, prepared = prepared)
}

staged_test_mesh = function(fixture, control, fixed = numeric()) {
  control$integration_mesh = prepare_staged_integration(fixture$theta,
    fixture$prepared, control, fixed)
  control
}

test_that('score refinement checks both normalizers and retains cellwise cancellation evidence', {
  # One positive cell: canonical mass1/score0, reference mass0.1/score10.
  # Canonical-normalized numerator disagreement is1, but the score gap is10.
  current = matrix(c(0,0), 1L, 2L)
  previous = matrix(c(10,-4), 1L, 2L)
  actual = staged_score_discrepancy(current, previous, 1, .1, c(0,0), c(10,-4))
  expect_equal(actual, c(10,4))
  expect_gt(actual[1L], abs(0-.1*10))
  # Opposing cell numerators cancel in the totals; their absolute evidence
  # remains in the gate even when both normalized total scores equalzero.
  current = matrix(c(100,-100), 2L, 1L)
  previous = matrix(0, 2L, 1L)
  expect_equal(staged_score_discrepancy(current, previous, c(.5,.5), c(.5,.5),
    0, 0), 100)
})

test_that('inner time budgets accumulate by person and exclude outer-only gaps', {
  tick = 0; now = function() tick
  prepared = list(group = c(1L,1L,2L), n_panels = 2L)
  budget = staged_inner_budget(prepared, list(seconds_per_panel = 5, seconds_total = 8), now)
  budget$begin(1L); tick = 2; budget$end()
  tick = 102
  expect_equal(budget$snapshot()$total_seconds, 2)
  budget$begin(2L); tick = 105; budget$end()
  budget$begin(3L); tick = 107; budget$end()
  expect_equal(budget$snapshot()$panel_seconds, c(5,2))
  expect_equal(budget$snapshot()$total_seconds, 7)
  tick = 500; budget$begin(2L); tick = 501
  error = tryCatch(budget$check(), panelselect_inner_time_error = identity)
  expect_s3_class(error, 'panelselect_inner_time_error')
  expect_match(conditionMessage(error), 'per-panel')
  expect_identical(error$diagnostics$panel, 1L)
  expect_equal(error$diagnostics$panel_seconds, c(6,2))
  budget$end(check_limit = FALSE)
  total = staged_inner_budget(prepared, list(seconds_per_panel = 10, seconds_total = 3), now)
  total$begin(1L); tick = 503; total$end()
  tick = 900; total$begin(3L); tick = 902
  expect_error(total$check(), 'total time budget', class = 'panelselect_inner_time_error')
  total$end(check_limit = FALSE)
  fresh = staged_inner_budget(prepared, list(seconds_per_panel = 5, seconds_total = 8), now)
  expect_equal(fresh$snapshot()$panel_seconds, c(0,0))
  expect_equal(fresh$snapshot()$total_seconds, 0)
})

test_that('active inner budgets bind direct, native GH and factored row mappings', {
  scenarios = list(c('probit_pln','gl'), c('pln_re','gh'), c('pln_re','agh'),
    c('probit_re_pln_re','gl'), c('probit_re_pln_re','gl'))
  fixtures = controls = fixed_values = vector('list', length(scenarios))
  for (i in seq_along(scenarios)) {
    scenario = scenarios[[i]]; fixture = pln_inner_fixture(scenario[1L])
    fixture$prepared$y = 1
    fixed = if (i == length(scenarios)) c(cor_re = 0, cor_error = 0) else numeric()
    fixture$theta[names(fixed)] = fixed
    control = if (scenario[2L] %in% c('gh','agh'))
      panel_control(quadrature = scenario[2L], nodes = 3L,
        observation_quadrature = 'gl', observation_nodes = 15L,
        observation_options = list(reference_nodes = 20L)) else
      staged_test_control('gl')
    control = staged_test_mesh(fixture, control, fixed)
    control$observation_options$seconds_per_panel = .001
    control$integration_mesh$binding = staged_binding(fixture$prepared, control, fixed)
    fixtures[[i]] = fixture; controls[[i]] = control; fixed_values[[i]] = fixed
  }
  original = staged_inner_budget
  local_mocked_bindings(staged_inner_budget = function(prepared, options) {
    tick = 0
    original(prepared, options, function() { tick <<- tick+1; tick })
  }, .package = 'PanelSelect')
  for (i in seq_along(scenarios)) {
    fixture = fixtures[[i]]; control = controls[[i]]; fixed = fixed_values[[i]]
    saved = serialize(control$integration_mesh, NULL)
    expect_error(evaluate_model(fixture$theta, fixture$prepared, control, fixed, score = FALSE),
      'per-panel time budget', class = 'panelselect_inner_time_error')
    expect_identical(serialize(control$integration_mesh, NULL), saved)
    expect_error(prepare_staged_integration(fixture$theta, fixture$prepared, control, fixed),
      'per-panel time budget', class = 'panelselect_inner_time_error')
  }
  # A factored outcome context retains original row IDs, including a skipped
  # unobserved row; both observed rows charge their original person's budget.
  fixture = core_fixture('selection_pln')
  fixture$prepared$.staged_rows = c(2L,4L,7L)
  fixed = c(cor_re = 0, cor_error = 0); fixture$theta[names(fixed)] = fixed
  control = staged_test_control('gl')
  control$integration_mesh = list(inner = vector('list', 7L))
  control$.staged_inner_budget = original(list(group = c(2L,1L,2L,1L,2L,2L,1L),
    n_panels = 2L), list(seconds_per_panel = 5, seconds_total = 8), function() 0)
  part = staged_panel_layout(fixture$theta, fixture$prepared, control, fixed)$parts$outcome
  context = make_model_context(fixture$theta, part, control, fixed)
  expect_identical(context$inner$row_ids, c(2L,7L))
  expect_identical(context$inner$budget, control$.staged_inner_budget)
})

test_that('nested native evaluation and analytic returns enforce an enclosing preparation clock', {
  fixture = pln_inner_fixture('pln_re'); fixture$prepared$y = 1
  control = staged_test_mesh(fixture, panel_control(quadrature = 'gh', nodes = 3L,
    observation_quadrature = 'gl', observation_nodes = 15L,
    observation_options = list(reference_nodes = 20L)))
  calls = 0L
  control$.staged_parent_clock = function() {
    calls <<- calls+1L
    integration_error('Enclosing staged preparation total time budget exhausted')
  }
  saved = serialize(control$integration_mesh, NULL)
  expect_error(evaluate_model(fixture$theta, fixture$prepared, control, score = FALSE),
    'Enclosing staged preparation', class = 'panelselect_integration_error')
  expect_identical(calls, 1L)
  expect_identical(serialize(control$integration_mesh, NULL), saved)
  # No inner callbacks or numerical outer nodes exist here. The successful
  # direct-kernel return must still invoke the enclosing overall guard.
  fixture = core_fixture('selection_probit'); fixture$prepared$structure = 'cross_sectional'
  fixture$prepared$model_id = 'probit_probit'
  fixture$theta = fixture$theta[!names(fixture$theta) %in% c('sd_re_outcome','sd_re_selection','cor_re')]
  fixture$prepared$group = 1:3; fixture$prepared$group_ptr = 0:3
  fixture$prepared$n_panels = 3L; fixture$prepared$id_values = as.character(1:3)
  control = staged_test_mesh(fixture, staged_test_control())
  control$.staged_parent_clock = function()
    integration_error('Enclosing analytic preparation total time budget exhausted')
  expect_error(evaluate_model(fixture$theta, fixture$prepared, control, score = FALSE),
    'Enclosing analytic preparation', class = 'panelselect_integration_error')
})

test_that('positive cell rules preserve volume through selective refinement and extension', {
  mesh = c(staged_cells(list(c(-2,0,2), c(-3,1,3))),
    list(method = 'gl', rule = staged_rule(15L, 20L), bound = 3))
  rule = staged_cell_rule(mesh)
  expect_equal(sum(exp(rule$log_weights)), 24, tolerance = 1e-12)
  refined = staged_refine_cells(mesh, 2L, 1L)
  expect_equal(sum(exp(staged_cell_rule(refined)$log_weights)), 24, tolerance = 1e-12)
  square = c(staged_cells(list(c(-3,0,3), c(-3,0,3))),
    list(method = 'gl', rule = staged_rule(15L, 20L), bound = 3))
  extended = staged_extend_cells(square, 5)
  expect_equal(sum(exp(staged_cell_rule(extended)$log_weights)), 100, tolerance = 1e-12)
})

test_that('bounded tensor templates preserve every cell rule', {
  for (q in 0:2) {
    cells = staged_cells(rep(list(c(-3,-.7,2)), q))
    plain = c(cells, list(method = 'gl', bound = 3,
      rule = staged_rule(12L, 15L)))
    cached = plain; cached$rule = staged_rule(12L, 15L, q)
    take = rev(seq_len(nrow(cells$lower)))
    for (reference in c(FALSE,TRUE)) {
      expected = staged_cell_rule(plain, take, reference)
      actual = staged_cell_rule(cached, take, reference)
      expect_identical(actual, expected)
      expect_identical(staged_cell_weights(cached, take, reference), expected$log_weights)
    }
  }
  higher = staged_rule(80L, 80L, 2L)
  expect_null(higher$tensor_template)
  expect_null(higher$reference$tensor_template)
})

test_that('outer batch allocation caps buffers while preserving complete rules', {
  fixture = core_fixture('selection_poisson')
  options = list(batch_nodes = 8192L)
  cells = staged_cells(rep(list(seq(-12,12,by = 1)), 2L))
  mesh = c(cells, list(method = 'gl', bound = 12,
    rule = staged_rule(15L, 20L, 2L)))
  ordinary = staged_outer_batch(mesh, fixture$prepared, options, 5L)
  expect_lte(max(ordinary$effective_nodes), options$batch_nodes)
  expect_lte(ordinary$estimated_buffer_bytes, ordinary$buffer_budget_bytes)
  expect_false(ordinary$one_cell_over_budget)
  expect_equal(ordinary$worst_panel_rows, 3L)
  wide = staged_outer_batch(mesh, fixture$prepared, options, 500L)
  expect_lte(wide$cells, ordinary$cells)
  expect_lte(wide$estimated_buffer_bytes, wide$buffer_budget_bytes)
  indivisible = staged_outer_batch(mesh, fixture$prepared, options, 100000L)
  expect_identical(indivisible$cells, 1L)
  expect_true(indivisible$one_cell_over_budget)
  expect_gt(indivisible$estimated_buffer_bytes, indivisible$buffer_budget_bytes)
  # A cell and each positive rule remain intact at every buffer setting.
  expect_equal(indivisible$effective_nodes[['canonical']], length(mesh$rule$nodes)^2)
  expect_equal(indivisible$effective_nodes[['reference']], length(mesh$rule[['reference']]$nodes)^2)
  expected_radius = max(abs(staged_cell_rule(mesh)$nodes),
    abs(staged_cell_rule(mesh, reference = TRUE)$nodes))
  expect_equal(staged_outer_radius(mesh), expected_radius, tolerance = 0)
  factored = list(route = 'independent', selection = mesh, outcome = mesh)
  expect_equal(staged_outer_radius(factored), staged_outer_radius(mesh), tolerance = 0)
  explicit = panel_control(quadrature = 'gl', observation_quadrature = 'gl',
    quadrature_options = list(batch_nodes = 256L), observation_options = list(batch_nodes = 64L))
  expect_equal(staged_options(explicit)$batch_nodes, 256L)
  expect_equal(staged_options(explicit, TRUE)$batch_nodes, 64L)
  gl = panel_control(quadrature = 'gl', observation_quadrature = 'gl')
  expect_equal(staged_options(gl)$batch_nodes, 8192L)
  expect_equal(staged_options(gl, TRUE)$batch_nodes, 128L)
})

test_that('batched inner diagnostics preserve scalar bounds and last-equal probes', {
  theta = c('selection.(Intercept)' = -.3, 'selection.x' = .1,
    'outcome.(Intercept)' = .25, 'outcome.x' = -.1, sd_re_selection = .6,
    sd_re_outcome = .7, cor_re = .35, sd_error = .5, cor_error = 0)
  mesh = list(binding = list(rows = 1L, panels = 1L,
    has_selection = TRUE, fixed = c(cor_error = 0), outcome_scale = 'sd_re_outcome'),
    anchor = theta)
  rule = list(theta = theta, row_ids = 1L, meshes = list(list(bound = 12)),
    projection = list(alpha = 1:2, beta = 3:4, x = matrix(c(1,-.2), 1L),
      w = matrix(c(1,.3), 1L), lambda_name = 'sd_re_outcome', lambda = .7,
      delta = .6, rho = .35, a = -.3, b = .25))
  B = c(.4,.8,.4,.5,.4,.8); A = c(.3,.5,.3,.4,.3,.5)
  y = rep(10, length(B)); log_factorial = lgamma(y+1)
  for (selected in c(FALSE,TRUE)) for (parameters in c(FALSE,TRUE))
    for (order in 0:1) {
      d = as.integer(selected)+1L+as.integer(parameters)+
        as.integer(selected && parameters)
      current = list(value = rep(log(.2), length(B)),
        gradient = outer(c(1,8,1,2,1,8), seq_len(d)/d))
      previous = list(value = log(c(.19,.199999,.19,.197,.19,.199999)),
        gradient = matrix(.1, length(B), d))
      scalar = staged_recorder(mesh, theta); scalar$outer_radius = 2
      for (i in seq_along(B)) staged_record_inner(scalar, rule, 1L,
        if (selected) A[i] else NULL, B[i], y[i], .5, if (selected) .25 else 0,
        log_factorial[i], list(value = current$value[i],
          gradient = if (order) current$gradient[i, ] else numeric()),
        list(value = previous$value[i],
          gradient = if (order) previous$gradient[i, ] else numeric()), order, parameters)
      for (batch in c(1L,2L,7L)) {
        actual = staged_recorder(mesh, theta); actual$outer_radius = 2
        staged_record_inner_batch(actual, rule, 1L,
          if (selected) A else NULL, B, y, .5, if (selected) .25 else 0,
          log_factorial, current, previous, order, parameters, batch_nodes = batch)
        for (field in c('inner_error','inner_score_error','inner_probes',
          'inner_score_probes','inner_score_probe_error'))
          expect_equal(actual[[field]], scalar[[field]], tolerance = 2e-14)
      }
      expect_equal(scalar$inner_probes[[1L]]$B, B[5L])
      if (order) expect_equal(scalar$inner_score_probes[[1L]]$B, B[6L])
      tails = vapply(seq_along(B), function(i) staged_inner_tail(
        if (selected) A[i] else NULL, B[i], y[i], .5,
        if (selected) .25 else 0, rule$meshes[[1L]], current$value[i],
        if (order) current$gradient[i, ] else numeric(d), parameters)$relative, numeric(1L))
      expect_equal(staged_inner_tail_relative_batch(y, rule$meshes[[1L]], current$value),
        tails, tolerance = 2e-14)
    }
  empty = staged_recorder(mesh, theta)
  staged_record_inner_batch(empty, rule, 1L, NULL, .4, 10, .5, 0,
    lgamma(11), list(value = -Inf), list(value = -Inf), 0L, FALSE)
  # Equal zero finite-rule masses have zero discrepancy; the omitted normal
  # tail remains a positive diagnostic. NA masses must retain their evidence.
  expect_true(is.finite(empty$inner_error[1L]))
  invalid = tryCatch(staged_record_inner_batch(staged_recorder(mesh, theta),
    rule, 1L, NULL, .4, 10, .5, 0, lgamma(11), list(value = NA_real_),
    list(value = -2), 0L, FALSE), panelselect_mesh_required = function(x) x)
  expect_s3_class(invalid, 'panelselect_mesh_required')
  expect_true(is.na(invalid$diagnostics$current_value))
  expect_equal(invalid$diagnostics$reference_value, -2)
})

test_that('staged GL likelihood agrees with the independently known probit beta anchor', {
  fixture = core_fixture('probit_re')
  fixture$theta[c('outcome.(Intercept)', 'outcome.x')] = 0
  fixture$theta['sd_re_outcome'] = 1
  # Under GL a probit RE panel takes moving cuts, so the staged engine is
  # called directly.
  control = staged_test_mesh(fixture, staged_test_control())
  value = evaluate_staged_integration(fixture$theta, fixture$prepared, control)
  k = sum(fixture$prepared$y); T = length(fixture$prepared$y)
  expect_equal(value$loglik, lbeta(k+1, T-k+1), tolerance = 1e-10)
  expect_equal(value$loglik, evaluate_staged_integration(fixture$theta, fixture$prepared,
    control, score = FALSE)$loglik, tolerance = 1e-12)
})

test_that('new count and coupled routes differentiate one immutable finite objective', {
  for (model in c('poisson_re', 'selection_poisson', 'selection_probit', 'selection_gaussian')) {
    fixture = core_fixture(model)
    control = staged_test_mesh(fixture, staged_test_control())
    saved = serialize(control$integration_mesh, NULL)
    actual = evaluate_model(fixture$theta, fixture$prepared, control)
    panel_errors = do.call(rbind, lapply(actual$diagnostics$panels,
      function(x) x$score_error))
    expect_true(all(panel_errors+1e-13 >= abs(actual$score-actual$reference_score)))
    # Disable acceptance gates only for the mathematical finite-objective
    # derivative test. The nodes/weights remain exactly those already checked.
    frozen = control; frozen$.staged_preparing = TRUE
    numerical = core_difference(function(theta) sum(evaluate_model(theta,
      fixture$prepared, frozen, score = FALSE)$loglik), fixture$theta)
    expect_lt(max(abs(colSums(actual$score)-numerical)), 5e-6)
    expect_identical(saved, serialize(control$integration_mesh, NULL))
    expect_equal(actual$loglik, evaluate_model(fixture$theta, fixture$prepared,
      control, score = FALSE)$loglik, tolerance = 1e-12)
    # Smaller outer batches regroup the same cells.
    small = control; small$quadrature_options$batch_nodes = 128L
    small$integration_mesh$binding = staged_binding(fixture$prepared, small, numeric())
    paired = evaluate_model(fixture$theta, fixture$prepared, small)
    expect_identical(paired$loglik, actual$loglik)
    expect_identical(paired$score, actual$score)
    expect_equal(paired$diagnostics$aggregate_loglik_error,
      actual$diagnostics$aggregate_loglik_error, tolerance = 1e-13)
    expect_equal(paired$diagnostics$aggregate_score_error,
      actual$diagnostics$aggregate_score_error, tolerance = 1e-13)
    replay = unserialize(serialize(control, NULL))
    expect_equal(evaluate_model(fixture$theta, fixture$prepared, replay)$score,
      actual$score, tolerance = 1e-12)
  }
})

test_that('PLN GL inner rules retain exact first through third finite derivatives', {
  fixture = pln_inner_fixture('probit_pln')
  theta = fixture$theta
  control = staged_test_mesh(fixture, staged_test_control('gl'))
  mesh = control$integration_mesh$inner[[1L]]
  point = c(a = theta[['selection.(Intercept)']], b = theta[['outcome.(Intercept)']],
    error = theta[['sd_error']], tau = theta[['cor_error']])
  value = function(x, order = 3L) staged_inner_value(x[['a']], x[['b']],
    fixture$prepared$y, x[['error']], x[['tau']], order,
    lgamma(fixture$prepared$y+1), mesh, probit_re_pln_re_inner_kernel, TRUE)
  actual = value(point)
  numerical = core_difference(function(x) value(x, 0L)$value, point)
  expect_lt(max(abs(drop(actual$gradient)-numerical)), 2e-6)
  for (j in 1:2) {
    plus = minus = point; step = 2e-5
    plus[j] = plus[j]+step; minus[j] = minus[j]-step
    d_first = (value(plus)$gradient-value(minus)$gradient)/(2*step)
    d_second = (value(plus)$hessian-value(minus)$hessian)/(2*step)
    expect_equal(drop(actual$hessian[1L, j, ]), drop(d_first), tolerance = 3e-6)
    expect_equal(matrix(actual$third[1L, , j, ], 2L), matrix(d_second, 2L), tolerance = 8e-6)
  }
})

test_that('paired PLN blocks reuse the evaluated reference without changing finite derivatives', {
  fixture = pln_inner_fixture('pln_re'); fixture$prepared$y = 1
  control = staged_test_mesh(fixture, staged_test_control('gl'))
  for (order in c(0L,1L,3L)) {
    create_rule = function() {
      copy = control; copy$.staged_recorder = staged_recorder(control$integration_mesh, fixture$theta)
      make_model_context(fixture$theta, fixture$prepared, copy, numeric())$inner
    }
    paired_rule = create_rule(); paired_rule$paired = TRUE
    canonical_rule = create_rule()
    reference_rule = create_rule(); reference_rule$reference = TRUE
    evaluate = function(rule) integrate_pln_inner_staged(NULL, c(-.3,.1,.9), rep(1,3L),
      fixture$theta[['sd_error']], 0, order, rep(0,3L), rule, pln_re_inner_kernel, TRUE)
    actual = evaluate(paired_rule); canonical = evaluate(canonical_rule); reference = evaluate(reference_rule)
    for (field in c('value','gradient','hessian','third')) {
      expect_equal(actual[[field]], canonical[[field]], tolerance = 1e-13)
      expect_equal(actual$reference_block[[field]], reference[[field]], tolerance = 1e-13)
    }
    expect_lt(sum(paired_rule$recorder$inner_work),
      sum(canonical_rule$recorder$inner_work)+sum(reference_rule$recorder$inner_work))
    expect_equal(paired_rule$recorder$inner_error, canonical_rule$recorder$inner_error, tolerance = 1e-13)
    expect_equal(paired_rule$recorder$inner_score_error, canonical_rule$recorder$inner_score_error, tolerance = 1e-13)
  }
})

test_that('tau-zero PLN error envelopes and acceptance are independent of outer batching', {
  fixture = core_fixture('selection_pln')
  fixed = c(cor_error = 0); fixture$theta[names(fixed)] = fixed
  large = staged_test_mesh(fixture, staged_test_control('gl'), fixed)
  saved = serialize(large$integration_mesh, NULL)
  actual = evaluate_model(fixture$theta, fixture$prepared, large, fixed)
  small = large; small$quadrature_options$batch_nodes = 128L
  small$integration_mesh$binding = staged_binding(fixture$prepared, small, fixed)
  paired = evaluate_model(fixture$theta, fixture$prepared, small, fixed)
  expect_identical(paired$loglik, actual$loglik)
  expect_identical(paired$score, actual$score)
  a = actual$diagnostics$panels[[1L]]; b = paired$diagnostics$panels[[1L]]
  for (field in c('outer_radius','nested_error','nested_score_error',
    'loglik_error','score_error','tail')) expect_equal(b[[field]], a[[field]], tolerance = 2e-12)
  expect_identical(staged_accuracy_reasons(b, staged_options(small), fixture$theta, fixed),
    staged_accuracy_reasons(a, staged_options(large), fixture$theta, fixed))
  expect_identical(saved, serialize(large$integration_mesh, NULL))
})

test_that('PLN outer adaptive GH uses row meshes without evaluation-history dependence', {
  fixture = pln_inner_fixture('pln_re')
  control = staged_test_budget(panel_control(quadrature = 'agh', nodes = 7L,
    observation_quadrature = 'gl', observation_nodes = 15L,
    observation_options = list(reference_nodes = 20L)))
  control = staged_test_mesh(fixture, control)
  saved = serialize(control$integration_mesh, NULL)
  frozen = control; frozen$.staged_preparing = TRUE
  first = evaluate_model(fixture$theta, fixture$prepared, frozen)
  probe = fixture$theta; probe['outcome.(Intercept)'] = probe['outcome.(Intercept)']+.1
  evaluate_model(probe, fixture$prepared, frozen)
  again = evaluate_model(fixture$theta, fixture$prepared, frozen)
  numerical = core_difference(function(theta) sum(evaluate_model(theta,
    fixture$prepared, frozen, score = FALSE)$loglik), fixture$theta)
  expect_equal(first$loglik, again$loglik, tolerance = 2e-12)
  expect_equal(first$score, again$score, tolerance = 2e-10)
  expect_match(first$diagnostics$nested_score_status, 'adaptive-geometry sensitivity')
  expect_match(first$diagnostics$reference_status, 'not enclosed')
  expect_lt(max(abs(colSums(first$score)-numerical)), 5e-6)
  expect_identical(saved, serialize(control$integration_mesh, NULL))
})

test_that('mesh binding includes count responses, offsets and fixed restrictions', {
  fixture = core_fixture('poisson_re')
  control = staged_test_mesh(fixture, staged_test_control())
  changed = fixture$prepared; changed$y[1L] = changed$y[1L]+1
  expect_error(evaluate_model(fixture$theta, changed, control), 'not bound')
  changed = fixture$prepared; changed$offset_outcome[1L] = .1
  expect_error(evaluate_model(fixture$theta, changed, control), 'not bound')
  expect_error(evaluate_model(fixture$theta, fixture$prepared, control,
    fixed = c(sd_re_outcome = .4)), 'not bound')
})

test_that('factored PLN inner error accounting retains the nonzero-rho pullback', {
  fixture = pln_inner_fixture('probit_re_pln_re')
  fixture$theta['cor_error'] = 0
  fixed = c(cor_error = 0)
  control = staged_test_mesh(fixture, staged_test_control('gl'), fixed)
  recorder = staged_recorder(control$integration_mesh, fixture$theta)
  recorder$outer_radius = 2
  control$.staged_recorder = recorder
  context = make_model_context(fixture$theta, fixture$prepared, control, fixed)
  current = list(value = log(.2), gradient = c(2,1))
  previous = list(value = log(.19), gradient = c(1.7,.8))
  staged_record_inner(recorder, context$inner, 1L, NULL, .4,
    fixture$prepared$y, fixture$theta[['sd_error']], 0,
    lgamma(fixture$prepared$y+1), current, previous, 1L, TRUE)
  rho = fixture$theta[['cor_re']]; lambda = fixture$theta[['sd_re_outcome']]
  expected = lambda*2*(1+abs(rho)/sqrt(1-rho^2))*abs(.2*2-.19*1.7)
  actual = exp(recorder$inner_score_error[1L, match('cor_re', names(fixture$theta))])
  expect_gte(actual+1e-15, expected)
  # The larger derivative discrepancy is at a different outer query from
  # the larger value discrepancy; both must reach between-stage refinement.
  staged_record_inner(recorder, context$inner, 1L, NULL, .6,
    fixture$prepared$y, fixture$theta[['sd_error']], 0,
    lgamma(fixture$prepared$y+1), current,
    list(value = log(.199999), gradient = c(.1,.1)), 1L, TRUE)
  probes = staged_inner_probes(recorder, 1L)
  expect_equal(length(probes), 2L)
  expect_equal(vapply(probes, function(x) x$B, numeric(1L)), c(.4,.6))
})

test_that('nested propagation uses actual finite GH moments and inactive inner rules skip staging', {
  fixture = pln_inner_fixture('pln_re')
  control = staged_test_budget(panel_control(quadrature = 'gh', nodes = 2L,
    observation_quadrature = 'gl', observation_nodes = 15L,
    observation_options = list(reference_nodes = 20L)))
  control = staged_test_mesh(fixture, control)
  actual = evaluate_model(fixture$theta, fixture$prepared, control)
  expect_equal(actual$diagnostics$finite_prior_mass, 1, tolerance = 2e-15)
  # GH2 has nodes +/-1 and E_H|Z|=1, exceeding sqrt(2/pi).
  expect_equal(actual$diagnostics$finite_prior_absolute_moment[1L,1L], 1, tolerance = 2e-15)
  expect_equal(actual$loglik, evaluate_model(fixture$theta, fixture$prepared,
    control, score = FALSE)$loglik, tolerance = 1e-12)
  fixture$theta['sd_error'] = 0
  expect_false(uses_staged_integration(fixture$prepared, control, c(sd_error = 0)))
})

test_that('analytic direct routes do not claim numerical quadrature checks', {
  fixture = core_fixture('selection_probit')
  fixture$theta = fixture$theta[!names(fixture$theta) %in% c('sd_re_outcome','sd_re_selection','cor_re')]
  fixture$prepared$structure = 'cross_sectional'; fixture$prepared$model_id = 'probit_probit'
  fixture$prepared$group = 1:3; fixture$prepared$group_ptr = 0:3
  fixture$prepared$n_panels = 3L; fixture$prepared$id_values = as.character(1:3)
  control = staged_test_mesh(fixture, staged_test_control())
  actual = evaluate_model(fixture$theta, fixture$prepared, control)
  expect_identical(actual$diagnostics$dimension, 0L)
  expect_false(actual$diagnostics$outer_order_checked)
  expect_false(actual$diagnostics$inner_order_checked)
  expect_false(actual$diagnostics$order_checked)
  expect_true(all(vapply(actual$diagnostics$panels,
    function(x) identical(x$nodes, 0L), logical(1L))))
  restricted = fixture; restricted$theta['cor_error'] = 0
  fixed = c(cor_error = 0)
  checked = staged_test_mesh(restricted, staged_test_control(), fixed)
  masked = evaluate_model(restricted$theta, restricted$prepared, checked, fixed)
  expect_identical(colnames(masked$score), names(restricted$theta))
  expect_equal(masked$score[, 'cor_error'], rep(0, restricted$prepared$n_panels))
})

test_that('staged probit likelihoods stay accurate when the posterior reaches the far tail', {
  # With index -3 + v, four successes put about 2.6e-6 of the first panel's
  # posterior mass beyond six SDs.
  data = data.frame(id = rep(1:2, each = 4), y = c(1, 1, 1, 1, 0, 1, 0, 1))
  prepared = prepare_data(y ~ 1, NULL, data, 'id', 'probit', 'fail')
  theta = c('outcome.(Intercept)' = -3, sd_re_outcome = 1)
  # Under GL a probit RE panel takes moving cuts, so the staged engine is
  # called directly.
  control = staged_test_mesh(list(theta = theta, prepared = prepared), staged_test_control())
  actual = evaluate_staged_integration(theta, prepared, control, score = FALSE)$loglik
  reference = function(y) {
    f = function(v) vapply(v, function(value) exp(sum(stats::pnorm(ifelse(y == 1, 1, -1)*(value - 3),
      log.p = TRUE)) + stats::dnorm(value, log = TRUE)), numeric(1))
    log(stats::integrate(f, -Inf, 3, rel.tol = 1e-12)$value + stats::integrate(f, 3, Inf, rel.tol = 1e-12)$value)
  }
  # The default per-panel log-likelihood tolerance is 1e-7.
  expect_lt(max(abs(actual - c(reference(c(1, 1, 1, 1)), reference(c(0, 1, 0, 1))))), 1e-7)
})
