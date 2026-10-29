score_tail_fixture = function(family = 'poisson', selected = FALSE) {
  one = matrix(1, 1L, 1L, dimnames = list(NULL, '(Intercept)'))
  theta = c('outcome.(Intercept)' = .25, sd_re_outcome = .7)
  if (selected) theta = c('selection.(Intercept)' = -.3, theta,
    sd_re_selection = .6, cor_re = .3)
  if (family == 'pln') theta = c(theta, sd_error = .5)
  if (selected && family == 'pln') theta = c(theta, cor_error = 0)
  prepared = list(model_id = if (selected) 'probit_re_pln_re' else
      if (family == 'pln') 'pln_re' else 'poisson_re',
    structure = 'panel', family = family, has_selection = selected,
    x = one, w = if (selected) one else matrix(numeric(), 1L, 0L),
    y = 1, d = 1, selected = TRUE, active = TRUE,
    offset_outcome = 0, offset_selection = 0, group = 1L,
    group_ptr = c(0L,1L), id_values = 'one', n_rows = 1L, n_panels = 1L)
  list(theta = theta[parameter_layout(prepared)$name], prepared = prepared)
}

# Staged GL person-effect rule with the observation rule `inner`. GL keeps its
# bound (12 by default), so an omitted tail it leaves unresolved stops.
score_tail_control = function(inner = 'gh') panel_control(quadrature = 'gl', nodes = 15L,
  observation_quadrature = inner, observation_nodes = 15L,
  quadrature_options = list(reference_nodes = 20L),
  observation_options = if (inner == 'gl') list(reference_nodes = 20L) else list())

# The interior discrepancy is exactly zero. Only a domain-dependent score
# tail rejects the initial bound; any interior subdivision is a wrong choice.
score_tail_outer_mock = function(calls, source = 'outer', coordinate = 'outcome.(Intercept)') {
  function(theta, prepared, control, fixed, mesh, recorder, clock, score = TRUE,
           outer_radius = NULL) {
    physical = identical(mesh$route, 'physical_pln')
    independent = identical(mesh$route, 'independent')
    bound = if (physical) min(mesh$u$bound, mesh$v$bound) else
      if (independent) min(mesh$selection$bound, mesh$outcome$bound) else mesh$bound
    inner_bound = if (length(control$integration_mesh$inner) &&
        !is.null(control$integration_mesh$inner[[1L]])) control$integration_mesh$inner[[1L]]$bound else Inf
    calls$outer = c(calls$outer, bound)
    if (source != 'outer') calls$inner = c(calls$inner, inner_bound)
    zero = setNames(numeric(length(theta)), names(theta)); tail_score = zero
    bad = if (source == 'outer') bound < 16 else inner_bound < 16
    if (bad && (source != 'combined_nested' || independent)) tail_score[coordinate] = 1
    nc = staged_interval_count(mesh)
    out = list(loglik = -1, score = zero, reference_score = zero,
      loglik_error = 0, score_error = tail_score, tail = 1e-12,
      outer_tail = 1e-12, inner_tail = 0,
      outer_tail_score = if (source == 'outer') tail_score else zero,
      inner_tail_score = if (source == 'inner') tail_score else zero,
      nested_error = 0,
      nested_score_error = if (source == 'combined_nested') tail_score else zero,
      cell_error = numeric(nc), cell_refinement_error = numeric(nc),
      dimension_error = matrix(0, nc, if (physical || independent) 2L else ncol(mesh$lower)),
      inner_probes = list(), nodes = staged_node_count(mesh), intervals = nc,
      dimension = if (physical || independent) 2L else ncol(mesh$lower))
    if (source == 'sum') {
      out$outer_tail_score = out$inner_tail_score = zero
      if (bound < 16) out$outer_tail_score[coordinate] = .6*5e-6
      if (inner_bound < 16) out$inner_tail_score[coordinate] = .6*5e-6
      out$score_error = out$outer_tail_score+out$inner_tail_score
    }
    if (independent) {
      out$inner_probes = list(list(row = 1L, A = NULL, B = .25, y = 1, log_factorial = 0))
      out$outcome = list(loglik = -1, inner_companion_log = c('1' = 0), inner_outer_mass = 1)
    }
    out
  }
}

score_tail_inner_mock = function(calls, amount = 1) {
  function(theta, prepared, control, fixed, mesh, probe, clock) {
    calls$probe = c(calls$probe, mesh$bound)
    d = if (is.null(probe$A)) 2L else 4L
    score = tail_score = numeric(d)
    if (mesh$bound < 16) tail_score[d] = amount
    list(finite = TRUE, loglik = 0, score = score,
      loglik_error = 0, score_error = tail_score, tail = 1e-15,
      tail_score = tail_score, cell_error = numeric(nrow(mesh$lower)))
  }
}

test_that('fixed outer coordinates do not stop the mesh and fixed-bound GL names score tails', {
  fixture = score_tail_fixture(); calls = new.env(parent = emptyenv()); calls$outer = numeric()
  local_mocked_bindings(staged_outer_panel = score_tail_outer_mock(calls),
    staged_refine_cells = function(...) stop('Unexpected interior refinement'), .package = 'PanelSelect')
  fixed = c('outcome.(Intercept)' = fixture$theta[['outcome.(Intercept)']])
  mesh = prepare_staged_integration(fixture$theta, fixture$prepared, score_tail_control(), fixed)
  expect_equal(calls$outer, 12)
  expect_equal(mesh$outer[[1L]]$bound, 12)
  calls$outer = numeric()
  error = tryCatch(prepare_staged_integration(fixture$theta, fixture$prepared,
    score_tail_control()), panelselect_mesh_required = identity)
  expect_s3_class(error, 'panelselect_mesh_required')
  expect_match(conditionMessage(error), 'omitted-score-tail.*maximum bound')
  expect_match(conditionMessage(error), 'raise quadrature_options$bound', fixed = TRUE)
  expect_equal(error$diagnostics$bound, 12)
  expect_gt(max(error$diagnostics$score_tail), max(error$diagnostics$score_budget))
  expect_equal(calls$outer, 12)
})

test_that('inner analytic score tail respects the original local and absolute budgets', {
  fixture = score_tail_fixture('pln'); probe = list(row = 1L, A = NULL,
    B = .25, y = 1, log_factorial = 0)
  calls = new.env(parent = emptyenv()); calls$probe = numeric()
  local_mocked_bindings(staged_inner_probe_value = score_tail_inner_mock(calls, 5e-6),
    staged_refine_cells = function(...) stop('Unexpected interior refinement'), .package = 'PanelSelect')
  control = score_tail_control('gl')
  seed_control = control; seed_control$.staged_seed = TRUE
  context = make_model_context(fixture$theta, fixture$prepared, seed_control, numeric())
  seed = staged_seed_mesh(context, fixture$prepared, fixture$theta, control, numeric(), inner = TRUE, row = 1L)
  saved = unserialize(serialize(seed, NULL))
  # The scaled score budget passes (5e-6<1e-5), so the row is accepted.
  accepted = staged_prepare_inner(fixture$theta, fixture$prepared, control, numeric(), seed,
    probe, function() NULL)
  expect_identical(accepted, seed)
  expect_equal(calls$probe, 12)
  # The existing absolute numerator target alone fails, while its probability
  # tail still passes; the fixed bound then stops the row.
  calls$probe = numeric()
  error = tryCatch(staged_prepare_inner(fixture$theta, fixture$prepared, control, numeric(),
    seed, probe, function() NULL, absolute_target = log(1e-7)), panelselect_mesh_required = identity)
  expect_s3_class(error, 'panelselect_mesh_required')
  expect_match(conditionMessage(error), 'omitted-score-tail.*maximum bound')
  expect_match(conditionMessage(error), 'raise observation_options$bound', fixed = TRUE)
  expect_equal(error$diagnostics$component, 'inner')
  expect_equal(calls$probe, 12)
  expect_identical(seed, saved)
  # The unchanged row gate requires all local derivatives, including error
  # when globally fixed; it is distinct from the panel's global free mask.
  fixed = c(sd_error = fixture$theta[['sd_error']])
  expect_error(staged_prepare_inner(fixture$theta, fixture$prepared, control, fixed, seed,
    probe, function() NULL, absolute_target = log(1e-7)), 'omitted-score-tail',
    class = 'panelselect_mesh_required')
})

test_that('an individually excessive fixed physical inner tail fails explicitly', {
  fixture = score_tail_fixture('pln', TRUE); fixed = c(cor_error = 0)
  calls = new.env(parent = emptyenv()); calls$outer = calls$inner = numeric()
  local_mocked_bindings(staged_outer_panel = score_tail_outer_mock(calls, 'inner', 'sd_error'),
    staged_refine_cells = function(...) stop('Unexpected interior refinement'), .package = 'PanelSelect')
  error = tryCatch(prepare_staged_integration(fixture$theta, fixture$prepared,
    score_tail_control('gl'), fixed), panelselect_mesh_required = identity)
  expect_s3_class(error, 'panelselect_mesh_required')
  expect_match(conditionMessage(error), 'omitted-inner-score-tail.*maximum bound')
  expect_match(conditionMessage(error), 'raise observation_options$bound', fixed = TRUE)
  expect_equal(error$diagnostics$bound, 12)
  expect_equal(calls$outer, 12)
})

test_that('an excessive physical outer score tail fails at the fixed GL bound', {
  fixture = score_tail_fixture('pln', TRUE); fixed = c(cor_error = 0)
  calls = new.env(parent = emptyenv()); calls$outer = calls$inner = numeric()
  local_mocked_bindings(staged_outer_panel = score_tail_outer_mock(calls, 'outer', 'sd_error'),
    staged_refine_cells = function(...) stop('Unexpected interior refinement'), .package = 'PanelSelect')
  error = tryCatch(prepare_staged_integration(fixture$theta, fixture$prepared,
    score_tail_control('gl'), fixed), panelselect_mesh_required = identity)
  expect_s3_class(error, 'panelselect_mesh_required')
  expect_match(conditionMessage(error), 'Physical PLN omitted-outer-score-tail.*maximum bound')
  expect_match(conditionMessage(error), 'raise quadrature_options$bound', fixed = TRUE)
  expect_equal(error$diagnostics$bound, c(u = 12, v = 12))
  expect_equal(calls$outer, 12)
})

test_that('two individually admissible physical score tails cannot pass with a shared budget exhausted', {
  fixture = score_tail_fixture('pln', TRUE); fixed = c(cor_error = 0)
  calls = new.env(parent = emptyenv()); calls$outer = calls$inner = numeric()
  local_mocked_bindings(staged_outer_panel = score_tail_outer_mock(calls, 'sum', 'sd_error'),
    staged_refine_cells = function(...) stop('Unexpected interior refinement'), .package = 'PanelSelect')
  # Each initial tail is 0.6 of the exact local score budget; their sum alone
  # rejects acceptance, and neither fixed domain can grow.
  error = tryCatch(prepare_staged_integration(fixture$theta, fixture$prepared,
    score_tail_control('gl'), fixed), panelselect_mesh_required = identity)
  expect_s3_class(error, 'panelselect_mesh_required')
  expect_match(conditionMessage(error), 'omitted-inner-score-tail.*maximum bound')
  expect_equal(calls$outer, 12)
})

test_that('independent combined nested tails recheck original outcome probes under absolute budgets and stop at the bound', {
  fixture = score_tail_fixture('pln', TRUE); fixture$theta['cor_re'] = 0
  fixed = c(cor_re = 0, cor_error = 0)
  calls = new.env(parent = emptyenv()); calls$outer = calls$inner = calls$probe = numeric()
  local_mocked_bindings(staged_outer_panel = score_tail_outer_mock(calls, 'combined_nested', 'sd_error'),
    staged_inner_probe_value = score_tail_inner_mock(calls, 5e-6),
    staged_refine_cells = function(...) stop('Unexpected interior refinement'), .package = 'PanelSelect')
  # The seed passes the scaled budget; the recheck of the recorded outcome
  # probe under its absolute budget fails at the fixed observation bound.
  error = tryCatch(prepare_staged_integration(fixture$theta, fixture$prepared,
    score_tail_control('gl'), fixed), panelselect_mesh_required = identity)
  expect_s3_class(error, 'panelselect_mesh_required')
  expect_match(conditionMessage(error), 'PLN observation omitted-score-tail.*maximum bound')
  expect_equal(error$diagnostics$component, 'inner')
  expect_equal(calls$probe, c(12,12))
})

test_that('real inner probes expose the unchanged analytic derivative-tail calculation', {
  fixture = score_tail_fixture('pln'); control = score_tail_control('gl')
  seed_control = control; seed_control$.staged_seed = TRUE
  context = make_model_context(fixture$theta, fixture$prepared, seed_control, numeric())
  mesh = staged_seed_mesh(context, fixture$prepared, fixture$theta, control, numeric(), inner = TRUE, row = 1L)
  probe = list(row = 1L, A = NULL, B = .25, y = 1, log_factorial = 0)
  value = staged_inner_probe_value(fixture$theta, fixture$prepared, control, numeric(), mesh,
    probe, function() NULL)
  maximum = exp(-1); derivative = exp(-1)+4*exp(-2)
  mass = maximum*2*pnorm(-mesh$bound)*exp(-value$loglik)
  expected = derivative*c(2*pnorm(-mesh$bound), 2*dnorm(mesh$bound))*
    exp(-value$loglik)+abs(value$score)*mass
  expect_equal(value$tail_score/expected, c(1,1), tolerance = 2e-12)
  expect_true(all(value$score_error >= value$tail_score))
})
