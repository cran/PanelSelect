test_that('GL and observation controls are explicit without changing GH defaults', {
  old = panel_control()
  expect_identical(old$quadrature, 'gh')
  expect_identical(old$nodes, 10L)
  expect_identical(old$observation_quadrature, 'gh')
  expect_identical(old$observation_nodes, 20L)
  outer = panel_control(quadrature = 'gl', quadrature_options = list(bound = 10))
  expect_identical(outer$nodes, 20L)
  expect_identical(outer$observation_nodes, 20L)
  expect_identical(validate_panel_control(outer), outer)
  inner = panel_control(observation_quadrature = 'gl')
  expect_identical(inner$nodes, 10L)
  expect_identical(inner$observation_nodes, 20L)
  gl = panel_control(observation_quadrature = 'gl', observation_nodes = 12,
    observation_options = list(reference_nodes = 15))
  expect_identical(gl$observation_options$reference_nodes, 15L)
  expect_error(panel_control(quadrature = 'gl', nodes = c(20, 20)), 'GL nodes')
  expect_error(panel_control(quadrature_options = list(bound = 8)), 'require')
  expect_error(panel_control(observation_options = list(bound = 8)), 'require')
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(max_cells = 1.5)), 'integer')
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(panel_score_tol = 0)), 'positive finite')
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(unknown = 1)), 'Unknown')
  retained = outer; retained$integration_mesh = list(epoch = 1L)
  retained$.rules = list(unused = TRUE)
  restored = validate_panel_control(retained)
  expect_identical(restored$integration_mesh, retained$integration_mesh)
  expect_null(restored$.rules)
})

# These doubles exercise restart/error ownership, rather than approximate an integral.
staged_interface_stage_fixture = function() {
  theta = c('outcome.(Intercept)' = 0)
  info = data.frame(name = names(theta), type = 'coefficient', fixed = FALSE)
  prepared = list(n_panels = 1L, n_rows = 2L, family = 'probit',
    model_id = 'probit_re', structure = 'panel')
  list(theta = theta, info = info, prepared = prepared, fixed = numeric(),
    control = fit_control(panel_control(quadrature = 'gl'), vcov = 'observed'))
}

staged_interface_stage_bindings = function(state, fail_phase = 'information', fail_until = 1L) {
  # Deliberate refinement requests arrive at the endpoint (the final evaluation,
  # which applies the stage's accuracy gate) or during information.
  refine = function(phase, theta) stop(staged_refinement_condition('Deliberate mesh refinement',
    diagnostics = list(test_reason = phase), probes = list(list(theta = theta))))
  list(
    prepare_staged_integration = function(theta, prepared, control, fixed,
        probes = NULL, previous = NULL, epoch = 1L) {
      state$epochs = c(state$epochs, epoch)
      state$controls[[epoch]] = control$quadrature_options
      if (epoch == 1L) expect_null(previous)
      list(epoch = epoch)
    },
    fit_evaluator = function(template, info, prepared, control, fixed) {
      best = NULL; checks = 0L
      evaluate = function(eta, trial = TRUE) {
        # The second non-trial evaluation of a stage is its endpoint.
        if (!trial) checks <<- checks + 1L
        if (!trial && checks == 2L && identical(fail_phase, 'endpoint_validation') &&
            control$integration_mesh$epoch <= fail_until) refine('endpoint_validation', eta)
        best <<- list(eta = eta, theta = eta)
        score = matrix(.2 - eta, 1L, 1L, dimnames = list(NULL, names(template)))
        list(theta = eta, loglik = -.5*(eta - .2)^2, score = score,
          optimizer_score = score, reference_score = score, valid = TRUE,
          diagnostics = list(method = 'test_immutable_positive_rule'))
      }
      list(evaluate = evaluate, rejected = function() 0L,
        rejections = function() list(total = 0L, by_class = integer(), by_reason = integer(), last = NULL),
        evaluations = function() checks, best_valid = function() best)
    },
    optimize_model = function(eta, evaluator, control, options) {
      state$optimizations = state$optimizations + 1L
      list(par = setNames(.2, names(eta)), counts = c('function' = 2L, gradient = 1L),
        convergence = 0L, message = NULL)
    },
    information_matrix = function(eta, evaluator, on_error = NULL) {
      if (identical(fail_phase, 'information') && tail(state$epochs, 1L) <= fail_until)
        refine('information', eta)
      diag(length(eta))
    }
  )
}

test_that('information refinement restarts optimization within the same three-stage budget', {
  fixture = staged_interface_stage_fixture()
  fixture$control$integration_mesh = list(epoch = 99L)
  state = new.env(parent = emptyenv())
  state$epochs = integer(); state$optimizations = 0L; state$controls = list()
  do.call(local_mocked_bindings, c(staged_interface_stage_bindings(state),
    list(.package = 'PanelSelect', .env = environment())))
  fit = do.call(fit_staged_model, fixture)
  expect_identical(state$epochs, 1:2)
  expect_identical(state$optimizations, 2L)
  expect_identical(fit$final$diagnostics$staging$stages_used, 2L)
  expect_identical(fit$final$diagnostics$staging$records[[1L]]$phase, 'information')
  expect_identical(fit$final$diagnostics$staging$records[[1L]]$error$diagnostics$test_reason, 'information')
  expect_identical(state$controls[[1L]], state$controls[[2L]])
  expect_identical(fit$control$integration_mesh$epoch, 2L)
  expect_true(fit$information_computed)
})

test_that('likelihood evaluations are recorded per stage and summed over stages', {
  fixture = staged_interface_stage_fixture()
  state = new.env(parent = emptyenv())
  state$epochs = integer(); state$optimizations = 0L; state$controls = list()
  do.call(local_mocked_bindings, c(staged_interface_stage_bindings(state,
    fail_phase = 'endpoint_validation', fail_until = 1L),
    list(.package = 'PanelSelect', .env = environment())))
  fit = do.call(fit_staged_model, fixture)
  records = fit$final$diagnostics$staging$records
  # Each stage's double evaluates its start and its endpoint (non-trial calls).
  expect_identical(vapply(records, function(x) x$evaluations, integer(1L)), c(2L, 2L))
  expect_identical(fit$evaluations, 4L)
})

test_that('unresolved meshes stop after three total stages and retain the causes', {
  fixture = staged_interface_stage_fixture()
  state = new.env(parent = emptyenv())
  state$epochs = integer(); state$optimizations = 0L; state$controls = list()
  do.call(local_mocked_bindings, c(staged_interface_stage_bindings(state,
    fail_phase = 'endpoint_validation', fail_until = 9L),
    list(.package = 'PanelSelect', .env = environment())))
  failure = tryCatch(do.call(fit_staged_model, fixture), error = identity)
  expect_s3_class(failure, 'panelselect_mesh_exhausted')
  expect_s3_class(failure, 'panelselect_mesh_required')
  expect_identical(state$epochs, 1:3)
  expect_identical(state$optimizations, 3L)
  expect_identical(failure$diagnostics$stages_used, 3L)
  expect_length(failure$diagnostics$stage_records, 3L)
  expect_identical(failure$diagnostics$cause$test_reason, 'endpoint_validation')
  for (record in failure$diagnostics$stage_records)
    expect_identical(record$error$diagnostics$test_reason, 'endpoint_validation')
})

test_that('the exhaustion message summarizes every stage and points to its records', {
  fixture = staged_interface_stage_fixture()
  state = new.env(parent = emptyenv())
  state$epochs = integer(); state$optimizations = 0L; state$controls = list()
  do.call(local_mocked_bindings, c(staged_interface_stage_bindings(state,
    fail_phase = 'endpoint_validation', fail_until = 9L),
    list(.package = 'PanelSelect', .env = environment())))
  failure = tryCatch(do.call(fit_staged_model, fixture), error = identity)
  message = conditionMessage(failure)
  expect_match(message, '^Integration mesh unresolved after three total stages: ')
  expect_match(message, paste0('By stage: stage 1 \\(endpoint_validation\\); ',
    'stage 2 \\(endpoint_validation\\); stage 3 \\(endpoint_validation\\)'))
  expect_match(message, '\\$diagnostics\\$stage_records', fixed = FALSE)
  records = list(list(epoch = 1L, phase = 'optimization',
    error = list(diagnostics = list(panel = 4L, reasons = c('likelihood refinement unresolved')))))
  expect_identical(staged_stage_summary(records),
    'stage 1 (optimization, panel 4: likelihood refinement unresolved)')
  # Accuracy-gate failures keep failing panels and per-panel reasons in probes.
  gate = list(list(epoch = 2L, phase = 'endpoint_validation', error = list(
    diagnostics = list(panels = list(list(loglik = 0))),
    probes = list(list(panels = c(2L, 5L), reasons = list(character(), 'tail unresolved',
      character(), character(), 'score refinement unresolved'))))))
  expect_identical(staged_stage_summary(gate),
    'stage 2 (endpoint_validation, panels 2, 5: tail unresolved; score refinement unresolved)')
})

test_that('endpoint and covariance failures consume one shared three-stage budget', {
  fixture = staged_interface_stage_fixture()
  state = new.env(parent = emptyenv())
  state$epochs = integer(); state$optimizations = 0L; state$controls = list()
  bindings = staged_interface_stage_bindings(state, fail_phase = 'endpoint_validation', fail_until = 1L)
  bindings$information_matrix = function(eta, evaluator, on_error = NULL) {
    if (tail(state$epochs, 1L) == 2L)
      stop(staged_refinement_condition('Covariance stencil needs refinement',
        diagnostics = list(stencil = eta), probes = list(list(theta = eta))))
    diag(length(eta))
  }
  do.call(local_mocked_bindings, c(bindings, list(.package = 'PanelSelect', .env = environment())))
  fit = do.call(fit_staged_model, fixture)
  expect_identical(state$epochs, 1:3)
  expect_identical(state$optimizations, 3L)
  expect_identical(vapply(fit$final$diagnostics$staging$records, function(x) x$phase, character(1L)),
    c('endpoint_validation', 'information', 'information'))
  expect_identical(fit$final$diagnostics$staging$stages_used, 3L)
  expect_match(fit$final$diagnostics$staging$records[[2L]]$error$message, 'Covariance stencil')
  expect_identical(state$controls[[1L]], state$controls[[3L]])
})

test_that('ordinary integration failures are preserved without mesh retries', {
  fixture = staged_interface_stage_fixture()
  calls = 0L
  local_mocked_bindings(prepare_staged_integration = function(...) {
    calls <<- calls + 1L
    integration_error('Deliberate resource exhaustion')
  }, .package = 'PanelSelect')
  failure = tryCatch(do.call(fit_staged_model, fixture), error = identity)
  expect_s3_class(failure, 'panelselect_integration_error')
  expect_match(conditionMessage(failure), 'resource exhaustion')
  expect_identical(calls, 1L)
  expect_identical(failure$diagnostics$stages_used, 1L)
})

test_that('observed-information stencils preserve requests for a new mesh', {
  eta = c('outcome.(Intercept)' = .2)
  evaluator = list(staged = TRUE, evaluate = function(point, trial = TRUE) {
    if (any(point != eta)) stop(staged_refinement_condition('Information stencil needs refinement',
      diagnostics = list(stencil = point), probes = list(list(theta = point))))
    list(loglik = 0, optimizer_score = matrix(0, 1L, 1L))
  }, information_steps = function(eta) 1e-3)
  failure = tryCatch(information_matrix(eta, evaluator), error = identity)
  expect_s3_class(failure, 'panelselect_mesh_required')
  expect_match(conditionMessage(failure), 'stencil needs refinement')
  expect_false(identical(failure$diagnostics$stencil, eta))
})

test_that('zero-dimensional models retain their reductions and serialize under GL', {
  data = data.frame(d = c(1, 0), y = c(1, NA))
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0, cor_error = .3)
  fit = panel_select(y ~ 1, d ~ 1, data, id = NULL, family = 'probit',
    fixed = theta, control = panel_control(quadrature = 'gl'))
  expect_identical(fit$model_id, 'probit_probit')
  expect_equal(fit$loglik_by_unit, c('1' = log(.25 + asin(.3)/(2*pi)), '2' = log(.5)), tolerance = 1e-12)
  expect_identical(model_dimension(theta, fit$prepared, fit$control, fit$fixed), 0L)
  expect_false(is.null(fit$integration$dimension))
  expect_true(all(fit$integration$dimension == 0L))
  restored = unserialize(serialize(fit, NULL))
  replay = evaluate_model(restored$coefficients, restored$prepared, restored$control, restored$fixed)
  # Saved contributions are named by observation; the evaluator's are not.
  expect_equal(replay$loglik, unname(fit$loglik_by_unit), tolerance = 1e-13)
  expect_equal(replay$score, `rownames<-`(fit$score, NULL), tolerance = 1e-13)
  updated_call = update(restored, evaluate = FALSE)
  expect_null(as.list(updated_call)$control$integration_mesh)
  updated = update(restored)
  expect_identical(updated$integration$staging$stages_used, 1L)
  expect_equal(predict(updated, data), predict(fit, data), tolerance = 1e-13)
  expect_error(predict(fit, data, conditioning = 'history', history = data), 'GL integration')
  expect_error(random_effects(fit), 'GL integration')
  expect_error(check_quadrature(fit), 'no numerical integral')
})

test_that('quadrature diagnostics use the dimensions of exact never-response reductions', {
  data = data.frame(id = rep(1:2, each = 3), d = 0, y = NA_real_)
  theta = c('selection.(Intercept)' = -.2, 'outcome.(Intercept)' = .3,
    sd_re_selection = 0, sd_re_outcome = .8, cor_re = 0)
  fit = panel_select(y ~ 1, d ~ 1, data, 'id', 'poisson',
    fixed = theta, control = panel_control(quadrature = 'gl'))
  expect_identical(model_dimension(theta, fit$prepared, fit$control, fit$fixed), 1L)
  expect_true(all(vapply(fit$integration$panels, function(x) x$dimension == 0L, logical(1L))))
  expect_equal(unname(fit$loglik_by_unit), rep(3*pnorm(.2, log.p = TRUE), 2L), tolerance = 1e-12)
  expect_error(check_quadrature(fit), 'no numerical integral')
})

test_that('serialized positive meshes replay and new-data updates prepare a fresh binding', {
  data = data.frame(id = rep(1:2, each = 3), y = c(0, 1, 1, 1, 0, 0))
  theta = c('outcome.(Intercept)' = 0, sd_re_outcome = 1)
  # Poisson RE panels are staged under GL; the reference integrates each
  # panel's likelihood over its person effect.
  reference = function(y) log(stats::integrate(function(v) vapply(v, function(value)
    prod(stats::dpois(y, exp(value))), numeric(1))*stats::dnorm(v), -Inf, Inf, rel.tol = 1e-12)$value)
  fit = panel_select(y ~ 1, NULL, data, 'id', 'poisson', fixed = theta,
    control = panel_control(quadrature = 'gl'))
  # The default per-panel log-likelihood tolerance is 1e-7.
  expect_lt(max(abs(fit$loglik_by_unit - c(reference(c(0, 1, 1)), reference(c(1, 0, 0))))), 1e-7)
  original_mesh = serialize(fit$control$integration_mesh, NULL)
  restored = unserialize(serialize(fit, NULL))
  replay = evaluate_model(restored$coefficients, restored$prepared, restored$control, restored$fixed)
  expect_equal(replay$loglik, unname(restored$loglik_by_unit), tolerance = 1e-12)
  updated = update(restored, data = data[1:2, ])
  expect_identical(updated$integration$staging$stages_used, 1L)
  expect_lt(abs(updated$loglik_by_unit - reference(c(0, 1))), 1e-7)
  expect_false(identical(updated$control$integration_mesh$binding, restored$control$integration_mesh$binding))
  expect_identical(serialize(restored$control$integration_mesh, NULL), original_mesh)
})
