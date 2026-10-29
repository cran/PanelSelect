# Staged-engine efficiency changes keep every finite rule and evaluation.

test_that('a hybrid evaluator builds its native caches once and refreshes transients', {
  local_test_seed(4)
  id = rep(1:6, each = 3L); x = stats::rnorm(18)
  y = stats::rpois(18, exp(.2 + .4*x + rep(stats::rnorm(6, sd = .5), each = 3L) +
    stats::rnorm(18, sd = .3)))
  prepared = prepare_data(y ~ x, NULL, data.frame(id = id, x = x, y = y), 'id', 'pln')
  info = parameter_layout(prepared); fixed = resolve_fixed(NULL, info, prepared)
  theta = c('outcome.(Intercept)' = .2, outcome.x = .4, sd_re_outcome = .5, sd_error = .3)
  control = staged_test_budget(panel_control(quadrature = 'agh', observation_quadrature = 'gl'))
  control$integration_mesh = prepare_staged_integration(theta, prepared, control, fixed)
  evaluator = fit_evaluator(theta, info, prepared, control, fixed)
  memo = environment(evaluator$evaluate)$prepared$integration_cache$memo
  eta = parameter_to_optimizer(theta, info)
  first = evaluator$evaluate(eta)
  built = lapply(c('staged_canonical', 'staged_reference'),
    function(key) memo[[key]]$integration_cache$memo)
  expect_true(all(vapply(built, is.environment, logical(1))))
  later = lapply(1:2, function(k) evaluator$evaluate(eta + .001*k))
  for (k in 1:2) expect_true(all(is.finite(later[[k]]$loglik)))
  # The same cache objects serve later evaluations: nothing was rebuilt.
  expect_identical(lapply(c('staged_canonical', 'staged_reference'),
    function(key) memo[[key]]$integration_cache$memo), built)
  # Each evaluation records its own observation work in a fresh recorder.
  expect_gt(first$diagnostics$inner_kernel_nodes, 0)
  expect_identical(first$diagnostics$inner_kernel_nodes, later[[1L]]$diagnostics$inner_kernel_nodes)
})

test_that('the staged workspace reuses mesh work without changing evaluations or fits', {
  local_test_seed(7)
  id = rep(1:8, each = 3L); x = stats::rnorm(24)
  y = stats::rpois(24, exp(.3 + .5*x + rep(stats::rnorm(8, sd = .6), each = 3L)))
  data = data.frame(id = id, x = x, y = y)
  prepared = prepare_data(y ~ x, NULL, data, 'id', 'poisson')
  info = parameter_layout(prepared); fixed = resolve_fixed(NULL, info, prepared)
  theta = c('outcome.(Intercept)' = .3, outcome.x = .5, sd_re_outcome = .6)
  control = staged_test_budget(panel_control(quadrature = 'gl'))
  control$integration_mesh = prepare_staged_integration(theta, prepared, control, fixed)
  evaluator = fit_evaluator(theta, info, prepared, control, fixed)
  workspace = environment(evaluator$evaluate)$prepared$staged_workspace
  expect_true(is.environment(workspace))
  for (shift in c(0, .002, -.003)) {
    eta = parameter_to_optimizer(theta, info) + shift
    warm = evaluator$evaluate(eta)
    direct = evaluate_model(parameter_from_optimizer(eta, theta, info), prepared, control, fixed)
    expect_identical(warm$loglik, direct$loglik)
    expect_identical(unname(warm$score), unname(direct$score))
  }
  expect_true(isTRUE(workspace$bound))
  expect_gt(workspace$bytes, 0)
  # Nothing environment-valued is stored with a staged fit or its mesh.
  fit = panel_select(y ~ x, NULL, data, 'id', 'poisson', start = theta, control = control)
  expect_null(fit$prepared$staged_workspace)
  has_environment = function(x) is.environment(x) ||
    (is.list(x) && any(vapply(x, has_environment, logical(1))))
  expect_false(has_environment(fit$control$integration_mesh))
})

test_that('meshes store one-dimensional rules and old template-carrying meshes still evaluate', {
  local_test_seed(7)
  id = rep(1:8, each = 3L); x = stats::rnorm(24)
  y = stats::rpois(24, exp(.3 + .5*x + rep(stats::rnorm(8, sd = .6), each = 3L)))
  prepared = prepare_data(y ~ x, NULL, data.frame(id = id, x = x, y = y), 'id', 'poisson')
  info = parameter_layout(prepared); fixed = resolve_fixed(NULL, info, prepared)
  theta = c('outcome.(Intercept)' = .3, outcome.x = .5, sd_re_outcome = .6)
  control = staged_test_budget(panel_control(quadrature = 'gl'))
  control$integration_mesh = prepare_staged_integration(theta, prepared, control, fixed)
  expect_true(all(vapply(control$integration_mesh$outer,
    function(mesh) is.null(mesh$rule$tensor_template), logical(1))))
  old = control
  old$integration_mesh$outer = lapply(old$integration_mesh$outer, function(mesh) {
    mesh$rule$tensor_template = staged_tensor_template(mesh$rule, ncol(mesh$lower)); mesh })
  expect_identical(evaluate_model(theta, prepared, old, fixed)[c('loglik', 'score')],
    evaluate_model(theta, prepared, control, fixed)[c('loglik', 'score')])
})

test_that('cross-sectional observation meshes keep no zero-dimension outer meshes', {
  local_test_seed(8)
  n = 30L; w = stats::rnorm(n); x = stats::rnorm(n)
  d = as.integer(.4 + .6*w + stats::rnorm(n) > 0)
  y = ifelse(d == 1, stats::rpois(n, exp(.2 + .3*x + stats::rnorm(n, sd = .4))), NA)
  data = data.frame(w = w, x = x, d = d, y = y)
  prepared = prepare_data(y ~ x, d ~ w, data, NULL, 'pln')
  info = parameter_layout(prepared)
  fixed = resolve_fixed(c(cor_error = 0), info, prepared)
  theta = c('selection.(Intercept)' = .4, selection.w = .6, 'outcome.(Intercept)' = .2,
    outcome.x = .3, sd_error = .4, cor_error = 0)[info$name]
  control = staged_test_budget(panel_control(quadrature = 'gl', observation_quadrature = 'gl'))
  control$integration_mesh = prepare_staged_integration(theta, prepared, control, fixed)
  expect_true(all(vapply(control$integration_mesh$outer, is.null, logical(1))))
  # A mesh saved with its zero-dimension outer meshes evaluates identically.
  old = control
  old$integration_mesh$outer = lapply(seq_len(prepared$n_panels), function(i)
    list(lower = matrix(numeric(), 1L, 0L), upper = matrix(numeric(), 1L, 0L), bound = 12))
  expect_identical(evaluate_model(theta, prepared, old, fixed)[c('loglik', 'score')],
    evaluate_model(theta, prepared, control, fixed)[c('loglik', 'score')])
})

test_that('saved staged fits keep per-panel summaries without per-cell vectors', {
  local_test_seed(7)
  id = rep(1:8, each = 3L); x = stats::rnorm(24)
  y = stats::rpois(24, exp(.3 + .5*x + rep(stats::rnorm(8, sd = .6), each = 3L)))
  fit = panel_select(y ~ x, NULL, data.frame(id = id, x = x, y = y), 'id', 'poisson',
    start = c('outcome.(Intercept)' = .3, outcome.x = .5, sd_re_outcome = .6),
    control = staged_test_budget(panel_control(quadrature = 'gl')))
  panel = fit$integration$panels[[1L]]
  expect_false(any(c('cell_error', 'cell_refinement_error', 'dimension_error') %in% names(panel)))
  expect_true(all(c('loglik', 'loglik_error', 'tail', 'nodes', 'intervals') %in% names(panel)))
})
