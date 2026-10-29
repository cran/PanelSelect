test_that('local observation partials differentiate the same scalar contributions', {
  rule = normal_rule(11L, 1L)
  examples = list(
    binary = list(point = c(.3, -.4, .25), r = 2L,
      fun = function(x) probit_re_probit_re_rows(x[1], x[2], 0, x[3], 3L)),
    gaussian = list(point = c(.3, -.4, .7, -.25), r = 2L,
      fun = function(x) linear_selection_log_contribution(x[1], x[2], .2, x[3], x[4], 3L, TRUE)),
    selected_pln = list(point = c(.3, -.4, .7, -.25), r = 2L,
      fun = function(x) integrate_pln_inner(x[1], x[2], 2, x[3], x[4], 3L,
        lgamma(3), rule, probit_re_pln_re_inner_kernel)),
    pln = list(point = c(-.4, .7), r = 1L,
      fun = function(x) integrate_pln_inner(NULL, x[1], 2, x[2], 0, 3L,
        lgamma(3), rule, pln_re_inner_kernel)))
  for (example in examples) {
    point = example$point; actual = example$fun(point)
    for (j in seq_along(point)) {
      plus = minus = point; h = 1e-5
      plus[j] = plus[j]+h; minus[j] = minus[j]-h
      upper = example$fun(plus); lower = example$fun(minus)
      expect_equal(actual$gradient[1, j], (upper$value-lower$value)/(2*h), tolerance = 2e-7)
      for (i in seq_len(example$r)) {
        expect_equal(actual$hessian[1, i, j],
          (upper$gradient[1, i]-lower$gradient[1, i])/(2*h), tolerance = 2e-7)
        for (k in seq_len(example$r))
          expect_equal(actual$third[1, i, k, j],
            (upper$hessian[1, i, k]-lower$hessian[1, i, k])/(2*h), tolerance = 2e-6)
      }
    }
  }
})

test_that('numeric panel kernels retain panel-fast node-major ordering', {
  for (model in core_models) {
    fixture = core_fixture(model); prepared = fixture$prepared
    index = c(1L, 2L, 3L, 1L, 3L)
    for (name in c('x', 'w')) prepared[[name]] = prepared[[name]][index, , drop = FALSE]
    for (name in c('y', 'd', 'offset_outcome', 'offset_selection', 'selected'))
      prepared[[name]] = prepared[[name]][index]
    prepared$group = c(1L, 1L, 1L, 2L, 2L)
    prepared$group_ptr = c(0L, 3L, 5L)
    prepared$n_rows = 5L; prepared$n_panels = 2L
    prepared$id_values = c('a', 'b')
    kernel = make_panel_kernel(fixture$theta, prepared, core_control(), numeric(),
                               order = 3L, latent = TRUE)
    z = matrix(seq_len(6L*kernel$q)/17-.3, 6L, kernel$q)
    actual = kernel$evaluate(z, prior = TRUE)
    for (latent in c(FALSE, TRUE)) {
      first = make_panel_kernel(fixture$theta, prepared, core_control(), numeric(),
        order = 1L, latent = latent)$evaluate(z, prior = TRUE)
      columns = if (latent) seq_len(ncol(actual$gradient)) else kernel$parameter_columns
      expect_equal(first$value, actual$value, tolerance = 1e-12)
      expect_equal(first$gradient, actual$gradient[, columns, drop = FALSE], tolerance = 1e-12)
    }
    effects = kernel$effects(z)
    for (block in 0:2) {
      rows = 1:2+2L*block
      expected = kernel$evaluate(z[rows, , drop = FALSE], prior = TRUE)
      expect_equal(actual$value[rows], expected$value, tolerance = 1e-12)
      expect_equal(actual$gradient[rows, , drop = FALSE], expected$gradient, tolerance = 1e-12)
      expect_equal(actual$hessian[rows, , , drop = FALSE], expected$hessian, tolerance = 1e-12)
      expect_equal(actual$third[rows, , , , drop = FALSE], expected$third, tolerance = 1e-12)
      single = kernel$effects(z[rows, , drop = FALSE])
      for (effect in c('u', 'v')) {
        expect_equal(effects[[effect]]$value[rows], single[[effect]]$value)
        expect_equal(effects[[effect]]$gradient[rows, , drop = FALSE], single[[effect]]$gradient)
      }
    }
  }
})

test_that('unused outcome rows do not enter local contributions or cached subsets', {
  fixture = core_fixture('selection_pln')
  prepared = prepare_integration(fixture$prepared, core_control())
  expected = make_panel_kernel(fixture$theta, prepared, core_control(), numeric())
  prepared$x[2, ] = NA_real_; prepared$y[2] = Inf
  prepared$model_data = prepared$integration_cache = NULL
  actual = make_panel_kernel(fixture$theta, prepared, core_control(), numeric())
  z = matrix(c(.1, -.2), 1L)
  expect_equal(actual$evaluate(z), expected$evaluate(z))
  prepared = prepare_integration(fixture$prepared, core_control())
  for (subset in list(subset_panels(prepared, 1L), selection_component(prepared))) {
    expect_null(subset$model_data)
    expect_null(subset$integration_cache)
  }
})

test_that('local normal derivatives agree with the reviewed tail fixture', {
  # Independent 1024-bit references; the former jet oracle repeated the
  # production closed form and its cancellation near -30.
  fixture = dget(testthat::test_path('fixtures', 'normal-tail.R'))$points
  take = match(c(-100, -40, -30-1e-6, -30+1e-6, 0, 8), fixture$x)
  expect_false(anyNA(take))
  actual = normal_log_derivatives(fixture$x[take], 3L)
  expect_equal(actual$value, fixture$value[take], tolerance = 1e-12)
  expect_equal(actual$first, fixture$first[take], tolerance = 1e-12)
  expect_equal(actual$second, fixture$second[take], tolerance = 1e-12)
  expect_equal(actual$third, fixture$third[take], tolerance = 1e-12)
})
