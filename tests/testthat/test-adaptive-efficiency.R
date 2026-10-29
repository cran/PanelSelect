# Compare the batched/reduced implementation with the frozen forward algorithm,
# not with regenerated expected results. Tight modes isolate algebra from stopping.
test_that('batched Newton directions retain scalar Cholesky safeguards', {
  for (q in 1:2) {
    H = array(0, c(6L, q, q)); H[, 1L, 1L] = c(1, 2, 0, -1, 1, 1)
    gradient = matrix(seq_len(6*q)/7, 6L, q); gradient[6L, ] = 0
    if (q == 2L) {
      H[, 1L, 2L] = c(.2, .3, 0, 0, 1, 0)
      H[, 2L, 1L] = H[, 1L, 2L]; H[2L, 2L, 1L] = -500
      H[, 2L, 2L] = c(1, 3, 2, 2, 1+.Machine$double.eps, 1)
    }
    expected = gradient
    for (i in 1:6) {
      factor = tryCatch(chol(matrix(H[i, , ], q, q)), error = function(e) NULL)
      step = if (is.null(factor)) gradient[i, ] else
        drop(backsolve(factor, forwardsolve(t(factor), gradient[i, ])))
      expected[i, ] = step/max(1, sqrt(sum(step^2))/5)
    }
    kernel = list(evaluate = function(z, prior)
      list(value = rep(0, 6L), gradient = gradient, hessian = -H))
    proposal = NULL
    values = list(evaluate = function(z, prior) {
      proposal <<- z
      stop('captured proposal')
    })
    expect_error(panel_mode(kernel, values, gradient*0,
      list(mode_maxit = 2L, mode_tol = 1e-10), 1e-12), 'captured proposal')
    expect_equal(proposal, expected, tolerance = 1e-12)
  }
})
test_that('batched geometry and reduced scores match the forward rule', {
  for (model in core_models) for (sign in c(-1, 1)) {
    fixture = core_fixture(model, sign)
    control = core_control('agh', nodes = 3L, inner = 7L)
    control$mode_tol = 1e-12
    prepared = prepare_integration(fixture$prepared, control)
    theta = fixture$theta
    old = reference_panel_geometry(theta, prepared, control, numeric())
    current = panel_geometry(theta, prepared, control, numeric())
    expect_equal(current[names(old)], old, tolerance = 1e-10, info = model)
    reference = with_mocked_bindings(
      reference_integrate_prepared(theta, prepared, control, numeric()),
      panel_geometry = reference_panel_geometry)
    actual = integrate_prepared(theta, prepared, control, numeric())
    expect_equal(actual$loglik, reference$loglik, tolerance = 1e-10, info = model)
    expect_equal(actual$score, reference$score, tolerance = 1e-10, info = model)
  }
})

test_that('reverse geometry contraction matches arbitrary forward differentials', {
  for (model in c('poisson_re', 'selection_probit', 'selection_gaussian')) {
    fixture = core_fixture(model)
    control = core_control('agh')
    prepared = prepare_integration(fixture$prepared, control)
    forward = panel_geometry(fixture$theta, prepared, control, numeric())
    reduced = panel_geometry(fixture$theta, prepared, control, numeric(), contracted = TRUE)
    q = ncol(forward$mode); p = length(fixture$theta)
    # Nonzero, asymmetric moments exercise the upper-triangle/off-diagonal rule.
    moments = matrix(seq_len(q + q*(q+1L)/2)/7 - .4, 1L)
    reference = forward$jacobian_gradient
    for (j in seq_len(p)) {
      reference[, j] = reference[, j] + sum(moments[, seq_len(q)]*forward$mode_gradient[, , j])
      index = q
      for (b in seq_len(q)) for (a in seq_len(b)) {
        index = index + 1L
        reference[, j] = reference[, j] + moments[, index]*forward$transform_gradient[, a, b, j]
      }
    }
    expect_equal(adaptive_score(reduced, moments), reference, tolerance = 1e-12)
  }
})

test_that('mode seeds are private, transactional and independent of trial order', {
  for (model in c('probit_re', 'poisson_re', 'selection_probit', 'selection_gaussian')) {
    fixture = core_fixture(model); info = parameter_layout(fixture$prepared)
    theta = fixture$theta; eta = parameter_to_optimizer(theta, info)
    control = panel_control(quadrature = 'agh', nodes = 3L)
    first = fit_evaluator(theta, info, fixture$prepared, control, numeric())
    other = fit_evaluator(theta, info, fixture$prepared, control, numeric())
    expected = first$evaluate(eta)
    changed = eta + .001*sin(seq_along(eta))
    a = first$evaluate(changed)
    b = other$evaluate(changed)
    expect_equal(a$loglik, b$loglik, tolerance = 1e-10)
    expect_equal(a$score, b$score, tolerance = 1e-10)
    expect_gt(first$mode_counts()[['warm']], 0L)
    expect_equal(other$mode_counts()[['warm']], 0L)
    returned = first$evaluate(eta)
    expect_equal(returned$score, expected$score, tolerance = 1e-10)
    expect_null(fixture$prepared$integration_cache$mode_state)
    modes = environment(first$evaluate)$modes
    seeds = modes$seeds
    original = evaluate_model
    with_mocked_bindings({
      rejected = first$evaluate(changed)
      expect_false(rejected$valid)
    }, evaluate_model = function(...) {
      result = original(...); result$score[] = NaN; result
    })
    expect_identical(modes$seeds, seeds)
    expect_identical(modes$pending, list())
    # A bad starting point falls back to the cold solve without committing it.
    for (key in names(modes$seeds)) modes$seeds[[key]]$mode[] = 1e300
    recovered = first$evaluate(eta)
    expect_true(recovered$valid)
    expect_equal(recovered$score, expected$score, tolerance = 1e-10)
    expect_gt(first$mode_counts()[['fallback']], 0L)
  }
})

test_that('information perturbations share and restore a fitted-point mode anchor', {
  fixture = core_fixture('selection_probit'); info = parameter_layout(fixture$prepared)
  theta = fixture$theta; eta = parameter_to_optimizer(theta, info)
  control = panel_control(quadrature = 'agh', nodes = 3L)
  evaluator = fit_evaluator(theta, info, fixture$prepared, control, numeric())
  evaluator$evaluate(eta)
  modes = environment(evaluator$evaluate)$modes
  seeds = modes$seeds
  information = information_matrix(eta, evaluator)
  expect_true(all(is.finite(information)))
  expect_null(modes$anchor)
  expect_identical(modes$seeds, seeds)
  cold_gradient = function(x) {
    value = parameter_from_optimizer(x, theta, info)
    score = evaluate_model(value, fixture$prepared, control)$score
    -colSums(score)*parameter_jacobian(value, info)
  }
  # The test compares warm and cold modes, so both use the evaluator's steps.
  cold = stats::optimHess(eta, function(x) 0, cold_gradient,
    control = list(ndeps = unname(evaluator$information_steps(eta))))
  expect_equal(information, cold, tolerance = 1e-8)
  with_mocked_bindings(
    expect_warning(information_matrix(eta, evaluator), 'Observed information is unavailable'),
    optimHess = function(...) stop('test failure'), .package = 'stats')
  expect_null(modes$anchor)
  expect_identical(modes$seeds, seeds)
})

test_that('distant optimizer trials do not reuse an unrelated mode', {
  fixture = core_fixture('probit_re'); info = parameter_layout(fixture$prepared)
  theta = fixture$theta; eta = parameter_to_optimizer(theta, info)
  evaluator = fit_evaluator(theta, info, fixture$prepared, panel_control(quadrature = 'agh', nodes = 3L), numeric())
  evaluator$evaluate(eta)
  nearby = eta; nearby[1L] = nearby[1L] + .001
  evaluator$evaluate(nearby)
  counts = evaluator$mode_counts()
  distant = eta; distant[1L] = distant[1L] + 4*(1 + abs(eta[1L]))
  result = evaluator$evaluate(distant)
  expect_true(result$valid)
  expect_equal(evaluator$mode_counts()[['warm']], counts[['warm']])
  expect_gt(evaluator$mode_counts()[['cold']], counts[['cold']])
})

test_that('PLN joint modes retain their cold-start branch policy', {
  fixture = core_fixture('selection_pln'); info = parameter_layout(fixture$prepared)
  theta = fixture$theta; eta = parameter_to_optimizer(theta, info)
  evaluator = fit_evaluator(theta, info, fixture$prepared,
    panel_control(quadrature = 'agh', nodes = 3L, observation_nodes = 7L), numeric())
  first = evaluator$evaluate(eta)
  evaluator$evaluate(eta + .001)
  second = evaluator$evaluate(eta)
  expect_identical(first$score, second$score)
  expect_equal(evaluator$mode_counts()[['warm']], 0L)
})

test_that('maxLik BFGS omits only its unused final Hessian', {
  skip_if_not_installed('maxLik')
  calls = list()
  local_mocked_bindings(maxLik = function(logLik, grad, start, method, control, ...) {
    calls[[method]] <<- list(...)
    list(estimate = start, code = 0L, iterations = 1L, message = 'test')
  }, .package = 'maxLik')
  evaluator = list(evaluate = function(x) list(loglik = -x^2, optimizer_score = matrix(-2*x, 1L)))
  for (method in c('BFGS', 'NR', 'BHHH'))
    optimize_model(.2, evaluator, fit_control(panel_control(optimizer = 'maxLik', method = method), vcov = 'observed'))
  expect_identical(calls$BFGS, list(finalHessian = FALSE))
  expect_identical(calls$NR, list())
  expect_identical(calls$BHHH, list())
})
