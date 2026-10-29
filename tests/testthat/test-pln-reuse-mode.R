test_that('normal reuse requires exact inputs and derivative orders', {
  original = normal_log_derivatives; calls = 0L
  local_mocked_bindings(normal_log_derivatives = function(x, order) {
    calls <<- calls+1L
    original(x, order)
  })
  cache = cache_normal_derivatives()
  x = c(-100, -30-1e-6, -30, -30+1e-6, 0, 8, Inf, NA_real_, NaN)
  for (order in 0:3) {
    expected = original(x, order)
    expect_identical(cache(x, order), expected)
    expect_identical(cache(x, order), expected)
  }
  expect_identical(calls, 4L)
  # A shared first element is only a prefilter, not a cache key.
  changed = x; changed[2L] = -20
  expect_identical(cache(changed, 1L), original(changed, 1L))
  expect_identical(calls, 5L)
  result = cache(x, 1L); result$first[1L] = 42
  changed[2L] = 10
  expect_identical(cache(x, 1L), original(x, 1L))
  expect_identical(cache(changed, 1L), original(changed, 1L))
  expect_identical(calls, 6L)
  # Signed zeros must not silently share a bitwise-different derivative result.
  for (value in list(c(-0, 1), c(0, 1), c(-0, 1)))
    expect_identical(cache(value, 3L), original(value, 3L))
  expect_identical(calls, 8L)
})

test_that('normal reuse obeys payload and entry limits without changing results', {
  original = normal_log_derivatives; calls = 0L
  local_mocked_bindings(normal_log_derivatives = function(x, order) {
    calls <<- calls+1L
    original(x, order)
  })
  # A scalar first-order entry holds input, value and first derivative: 24 bytes.
  cache = cache_normal_derivatives(limit = 24, max_entries = 10L)
  for (value in c(1, 1, 2, 2, 1))
    expect_identical(cache(value, 1L), original(value, 1L))
  expect_identical(calls, 3L)
  cache = cache_normal_derivatives(limit = 1024, max_entries = 1L)
  calls = 0L
  for (value in c(1, 2, 1, 2))
    expect_identical(cache(value, 1L), original(value, 1L))
  expect_identical(calls, 3L)
  cache = cache_normal_derivatives(limit = 0)
  calls = 0L
  for (i in 1:2) expect_identical(cache(1, 1L), original(1, 1L))
  expect_identical(calls, 2L)
})

test_that('PLN inner mode partials retain the full index derivatives', {
  A = c(-4, .3, 5); B = c(-.3, .7, 2); y = c(0, 4, 25)
  rule = normal_rule(9L, 1L)
  cached_normal = cache_normal_derivatives()
  cached_kernel = function(...) probit_re_pln_re_inner_kernel(...,
    normal_evaluate = cached_normal)
  kernels = list(pln_re_inner_kernel, probit_pln_inner_kernel,
    probit_re_pln_re_inner_kernel, cached_kernel)
  for (i in seq_along(kernels)) for (tau in c(-.4, 0)) for (order in 0:3) {
    arguments = list(A = if (i == 1L) NULL else A, B = B, y = y,
      gamma = .8, tau = tau, order = order, log_factorial = lgamma(y+1),
      rule = rule, kernel = kernels[[i]])
    full = do.call(integrate_pln_inner, arguments)
    mode = do.call(integrate_pln_inner, c(arguments, list(parameters = FALSE)))
    if (i == 4L) {
      uncached = arguments; uncached$kernel = probit_re_pln_re_inner_kernel
      expect_identical(full, do.call(integrate_pln_inner, uncached))
    }
    # The free tau score survives at tau=0; only a mode request omits its column.
    if (i > 1L && tau == 0 && order)
      expect_true(any(abs(full$gradient[, 4L]) > 1e-5))
    take = seq_len(full$r)
    full$coordinates = full$coordinates[take]
    if (order) full$gradient = full$gradient[, take, drop = FALSE]
    if (order >= 2L) full$hessian = full$hessian[, , take, drop = FALSE]
    if (order >= 3L) full$third = full$third[, , , take, drop = FALSE]
    expect_identical(mode, full)
  }
})

test_that('PLN panel mode requests omit only parameter derivative columns', {
  cases = list(
    list(model = 'pln_re', fixed = numeric()),
    list(model = 'pln_re', fixed = c(sd_error = 0)),
    list(model = 'selection_pln', fixed = numeric()),
    list(model = 'selection_pln', fixed = numeric(), tau = 0),
    list(model = 'selection_pln', fixed = c(cor_error = 0)),
    list(model = 'selection_pln', fixed = c(sd_error = 0, cor_error = 0)))
  control = core_control('agh', nodes = 3L, inner = 7L)
  for (case in cases) {
    fixture = core_fixture(case$model)
    fixture$theta[names(case$fixed)] = case$fixed
    if (!is.null(case$tau)) fixture$theta['cor_error'] = case$tau
    prepared = prepare_integration(fixture$prepared, control)
    q = model_dimension(fixture$theta, prepared, control, case$fixed)
    z = matrix(sin(seq_len(3L*q)), 3L, q)
    for (order in 0:3) {
      full = make_panel_kernel(fixture$theta, prepared, control, case$fixed,
        order, parameters = TRUE, latent = TRUE)$evaluate(z, prior = TRUE)
      mode = make_panel_kernel(fixture$theta, prepared, control, case$fixed,
        order, parameters = FALSE, latent = TRUE)$evaluate(z, prior = TRUE)
      take = seq_len(q)
      full$gradient = full$gradient[, take, drop = FALSE]
      if (order >= 2L) full$hessian = full$hessian[, , take, drop = FALSE]
      if (order >= 3L) full$third = full$third[, , , take, drop = FALSE]
      expect_equal(mode, full, tolerance = 1e-12, info = case$model)
    }
  }
})

test_that('fixed PLN normal reuse remains local to each parameter evaluation', {
  fixture = core_fixture('selection_pln')
  control = core_control(nodes = c(3L, 5L), inner = 7L)
  prepared = prepare_integration(fixture$prepared, control)
  original = cache_normal_derivatives; caches = 0L
  local_mocked_bindings(cache_normal_derivatives = function(...) {
    caches <<- caches+1L
    original(...)
  })
  changed = fixture$theta; changed['cor_error'] = 0
  points = list(fixture$theta, changed, fixture$theta)
  before = serialize(list(fixture, control), NULL)
  actual = lapply(points, function(theta) evaluate_model(theta, prepared, control))
  expect_identical(caches, 3L)
  expected = with_mocked_bindings(
    lapply(points, function(theta) evaluate_model(theta, prepared, control)),
    cache_normal_derivatives = function(...) normal_log_derivatives)
  expect_identical(actual, expected)
  expect_identical(actual[[1L]], actual[[3L]])
  expect_identical(serialize(list(fixture, control), NULL), before)
  expect_null(prepared$selection_normal)
  expect_null(prepared$integration_cache$selection_normal)
  # Fixed zero and undifferentiated evaluations do not allocate this cache.
  evaluate_model(changed, prepared, control, fixed = c(cor_error = 0))
  evaluate_model(fixture$theta, prepared, control, score = FALSE)
  expect_identical(caches, 3L)
})
