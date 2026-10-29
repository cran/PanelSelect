# RE is the default start: family-specific auxiliaries for panel selection and
# silent skips for fixed zero scales and unidentified auxiliaries. Fits stored before init
# existed keep pooled starts (test-interface-inputs.R).

default_prepared = function(seed = 2301, panels = 20L, waves = 4L, family = 'poisson') {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + stats::rnorm(n) > 0)
  data$y = stats::rpois(n, exp(.2 + .4*data$x + .5*u[data$id]))
  if (family == 'probit') data$y = as.integer(data$y > 1)
  data$y[data$d == 0] = NA
  prepared = prepare_data(y ~ x, d ~ w, data, 'id', family)
  info = parameter_layout(prepared)
  list(data = data, prepared = prepared, info = info)
}

starts = function(f, fixed = numeric(), control = fit_control()) {
  info = f$info
  fixed = resolve_fixed(fixed, info, f$prepared)
  info$fixed = info$name %in% names(fixed)
  initial_parameters(f$prepared, info, NULL, fixed, control)
}

test_that('RE is the default start and equals the explicit re path', {
  expect_identical(eval(formals(panel_select)$init)[1L], 'default')
  expect_identical(fit_control()$init, 'default')
  f = default_prepared()
  default = starts(f)
  expect_identical(default, starts(f, control = fit_control(init = 're')))
  expect_identical(default, starts(f, control = fit_control(init = 'default')))
  expect_false(identical(default, starts(f, control = fit_control(init = 'pooled'))))
})

test_that('initialization aliases preserve starts, fitted models and updates', {
  f = default_prepared()
  control = panel_control(nodes = 4L, maxit = 0L)
  fit = function(method) panel_select(y ~ x, d ~ w, f$data, 'id', 'poisson',
    control = control, init = method)
  reference = list(default = fit('default'), pooled = fit('pooled'))
  for (alias in c('re', 'cross_sectional')) {
    canonical = if (alias == 're') 'default' else 'pooled'
    actual = fit(alias)
    expected = reference[[canonical]]
    expect_identical(actual$model_id, 'probit_re_poisson_re')
    expect_identical(actual$control$init, canonical)
    expect_identical(coef(actual), coef(expected))
    expect_identical(actual$loglik_by_unit, expected$loglik_by_unit)
    expect_identical(actual$score, expected$score)
    expect_identical(update(actual, init = alias, evaluate = FALSE)$init, canonical)
    restarted = update(actual, init = alias, start = NULL)
    expect_identical(coef(restarted), coef(expected))
    # Saved spellings and carried controls must not change initialization or
    # cause a warning when they name the same method as the public argument.
    saved = unserialize(serialize(actual, NULL))
    saved$control$init = alias
    expect_silent({ updated = update(saved, start = NULL) })
    expect_identical(coef(updated), coef(expected))
    expect_identical(updated$control$init, canonical)
    expect_silent({ carried = panel_select(y ~ x, d ~ w, f$data, 'id', 'poisson',
      control = saved$control, init = canonical) })
    expect_identical(coef(carried), coef(expected))
  }
})

test_that('fixed zero scales and one-row probit panels skip auxiliaries without a warning', {
  calls = 0L
  local_mocked_bindings(panel_select = function(...) {
    calls <<- calls + 1L
    stop('auxiliary called')
  }, .package = 'PanelSelect')
  f = default_prepared()
  zero_scales = c(sd_re_selection = 0, sd_re_outcome = 0)
  expect_silent({ zero = starts(f, zero_scales) })
  expect_identical(zero, starts(f, zero_scales, fit_control(init = 'pooled')))
  single = default_prepared(panels = 40L, waves = 1L, family = 'probit')
  expect_silent({ one = starts(single, c(cor_re = .5)) })
  expect_identical(one, starts(single, c(cor_re = .5), fit_control(init = 'pooled')))
  expect_identical(calls, 0L)
})
