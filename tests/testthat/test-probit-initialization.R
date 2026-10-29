probit_initial_fixture = function(family = "probit", cross = FALSE, selected = TRUE,
                                  no_intercept = FALSE) {
  local_test_seed(2402)
  n = 20L; waves = 4L
  data = data.frame(id = rep(seq_len(n), each = waves), x = rnorm(n*waves),
    w = rnorm(n*waves), group = factor(rep(c("a", "b"), n*waves/2)),
    off = seq(-.2, .2, length.out = n*waves))
  u = rep(rnorm(n), each = waves); v = rep(rnorm(n), each = waves)
  data$d = as.numeric(.2 + .4*data$w + .8*u + rnorm(n*waves) > 0)
  data$y = as.numeric(-.1 + .3*data$x + .7*v + rnorm(n*waves) > 0)
  if (family == "gaussian") data$y = data$x + v + rnorm(n*waves)
  if (family %in% c("poisson", "pln")) data$y = rpois(n*waves, exp(.2 + .2*data$x + .3*v))
  if (selected) data$y[data$d == 0] = NA_real_
  outcome = if (no_intercept) y ~ 0 + offset(off) else y ~ poly(x, 2) + group + offset(off)
  selection = if (!selected) NULL else if (no_intercept) d ~ 0 + offset(off) else
    d ~ poly(w, 2) + group + offset(off)
  prepared = prepare_data(outcome, selection, data, if (cross) NULL else "id", family)
  list(data = data, formula = outcome, selection = selection, prepared = prepared,
    info = parameter_layout(prepared))
}

test_that('GL initialization propagates numerical controls and creates fresh auxiliary stages', {
  f = probit_initial_fixture()
  control = fit_control(panel_control(quadrature = 'gl', maxit = 0,
    quadrature_options = list(bound = 10, max_nodes = 2000000)), init = 're')
  control$integration_mesh = list(epoch = 99L)
  controls = list()
  local_mocked_bindings(panel_select = function(outcome, selection, data, id, family,
      start = NULL, fixed = NULL, control = panel_control(), subset = NULL, na_action = 'fail',
      init = 're', vcov = 'opg') {
    prepared = prepare_data(outcome, selection, data, id, family)
    info = parameter_layout(prepared)
    theta = setNames(rep(.3, nrow(info)), info$name)
    theta[info$type == 'sd'] = .8
    controls[[length(controls) + 1L]] <<- c(control, list(.init = init, .vcov = vcov))
    list(coefficients = theta, optimizer = list(converged = TRUE))
  }, .package = 'PanelSelect')
  initial_parameters(f$prepared, f$info, NULL, numeric(), control)
  expect_length(controls, 2L)
  for (aux in controls) {
    expect_identical(aux$quadrature, 'gl')
    expect_identical(aux$nodes, 20L)
    expect_identical(aux$quadrature_options, control$quadrature_options)
    expect_null(aux$integration_mesh)
    expect_identical(aux$optimizer, 'optim')
    expect_identical(aux$.init, 'pooled')
    expect_identical(aux$.vcov, 'opg')
    expect_identical(aux$maxit, panel_control()$maxit)
  }
})

test_that('GL point-cap defaults apply at both levels without changing explicit budgets', {
  gl = panel_control(quadrature = 'gl', observation_quadrature = 'gl')
  for (options in list(gl$quadrature_options, gl$observation_options)) {
    expect_identical(options$max_nodes, 100000000L)
    expect_identical(options$seconds_per_panel, 60)
    expect_identical(options$seconds_total, 300)
    expect_identical(options$panel_loglik_tol, 1e-7)
    expect_identical(options$panel_score_tol, 1e-5)
    expect_identical(options$tail_tol, 1e-9)
  }
  explicit_gl = panel_control(quadrature = 'gl', quadrature_options = list(max_nodes = 300),
    observation_quadrature = 'gl', observation_options = list(max_nodes = 700))
  expect_identical(explicit_gl$quadrature_options$max_nodes, 300L)
  expect_identical(explicit_gl$observation_options$max_nodes, 700L)
  expect_identical(panel_control()$quadrature, 'gh')
  expect_identical(panel_control()$observation_quadrature, 'gh')
})

test_that('modern binary auxiliary nonconvergence preserves only the failed equation starts', {
  f = probit_initial_fixture()
  supplied = c('selection.poly(w, 2)1' = .35, sd_re_outcome = .62)
  pooled = initial_parameters(f$prepared, f$info, supplied, numeric(), fit_control(init = 'pooled'))
  for (failed in c('selection', 'outcome')) {
    calls = 0L
    local_mocked_bindings(panel_select = function(outcome, selection, data, id, family, ...) {
      calls <<- calls + 1L
      info = parameter_layout(prepare_data(outcome, NULL, data, id, family))
      value = setNames(rep(.47, nrow(info)), info$name)
      list(coefficients = value, optimizer = list(converged =
        calls != match(failed, c('selection', 'outcome'))))
    }, .package = 'PanelSelect')
    expect_warning({actual = initial_parameters(f$prepared, f$info, supplied, numeric(),
      fit_control(init = 're'))}, paste(failed, 'RE initialization failed'))
    failed_names = if (failed == 'selection') c(paste0('selection.', colnames(f$prepared$w)),
      'sd_re_selection') else c(paste0('outcome.', colnames(f$prepared$x)), 'sd_re_outcome')
    expect_identical(actual[failed_names], pooled[failed_names])
    expect_identical(actual[names(supplied)], supplied)
    successful = setdiff(c(if (failed == 'selection') paste0('outcome.', colnames(f$prepared$x)) else
      paste0('selection.', colnames(f$prepared$w)), if (failed == 'selection') 'sd_re_outcome' else
      'sd_re_selection'), names(supplied))
    expect_equal(unname(actual[successful]), rep(.47, length(successful)))
  }
})

test_that("binary RE starts reuse prepared matrices and independent auxiliary controls", {
  fixture = probit_initial_fixture()
  data = fixture$data; data$id[c(2L, 7L)] = NA_integer_
  prepared = prepare_data(fixture$formula, fixture$selection, data, "id", "probit", "omit")
  info = parameter_layout(prepared)
  start = c("selection.poly(w, 2)1" = .8, sd_re_selection = .55)
  fixed = c("selection.(Intercept)" = .4, sd_re_outcome = .9)
  info$fixed = info$name %in% names(fixed)
  control = fit_control(panel_control(quadrature = "gh", nodes = c(3, 5),
    optimizer = "maxLik", method = "NR", maxit = 0, reltol = 1e-4,
    mode_maxit = 17, mode_tol = 1e-7), init = "re", vcov = "sandwich")
  pooled_control = control; pooled_control$init = "pooled"
  pooled = initial_parameters(prepared, info, start, fixed, pooled_control)
  calls = list()
  local_mocked_bindings(panel_select = function(outcome, selection, data, id, family,
      start = NULL, fixed = NULL, control = panel_control(), subset = NULL, na_action = "fail",
      init = "re", vcov = "opg") {
    p = prepare_data(outcome, selection, data, id, family)
    layout = parameter_layout(p)
    value = setNames(seq_len(nrow(layout))/10 + length(calls), layout$name)
    value[layout$type == "sd"] = .6 + length(calls)/10
    value[names(fixed)] = fixed
    calls[[length(calls) + 1L]] <<- list(prepared = p, start = start, fixed = fixed,
      control = control, init = init, vcov = vcov, coefficients = value)
    list(coefficients = value, optimizer = list(converged = TRUE))
  }, .package = "PanelSelect")
  result = initial_parameters(prepared, info, start, fixed, control)
  expect_length(calls, 2L)
  designs = list(prepared$w, prepared$x[prepared$selected, , drop = FALSE])
  responses = list(prepared$d, prepared$y[prepared$selected])
  offsets = list(prepared$offset_selection, prepared$offset_outcome[prepared$selected])
  blocks = list(c(paste0("selection.", colnames(prepared$w)), "sd_re_selection"),
    c(paste0("outcome.", colnames(prepared$x)), "sd_re_outcome"))
  for (j in 1:2) {
    call = calls[[j]]
    expect_equal(unname(call$prepared$x), unname(designs[[j]]), tolerance = 0)
    expect_equal(call$prepared$y, responses[[j]])
    expect_equal(call$prepared$offset_outcome, offsets[[j]])
    expect_false(call$prepared$has_selection)
    expect_identical(call$prepared$family, "probit")
    expect_identical(call$control$quadrature, "gh")
    expect_equal(call$control$nodes, c(3, 5)[j])
    expect_identical(call$control$mode_maxit, 17L)
    expect_identical(call$control$mode_tol, 1e-7)
    expect_identical(call$control$optimizer, "optim")
    expect_identical(call$control$method, "BFGS")
    expect_identical(call$vcov, "opg")
    expect_identical(call$init, "pooled")
    expect_false(any(c("na_action", "init", "vcov") %in% names(call$control)))
    expect_identical(call$control$maxit, panel_control()$maxit)
    expect_identical(call$control$reltol, panel_control()$reltol)
    expect_equal(unname(call$start), unname(pooled[blocks[[j]]]))
    expected = unname(call$coefficients)
    names(expected) = blocks[[j]]
    expected[intersect(names(expected), names(fixed))] = fixed[intersect(names(expected), names(fixed))]
    expected[intersect(names(expected), names(start))] = start[intersect(names(expected), names(start))]
    expect_equal(result[blocks[[j]]], expected)
    expect_equal(unname(call$fixed), unname(fixed[intersect(blocks[[j]], names(fixed))]))
  }
  expect_identical(result[names(start)], start)
  expect_identical(result[names(fixed)], fixed)
  expect_equal(result[c("cor_re", "cor_error")], pooled[c("cor_re", "cor_error")])
})

test_that("complete starts and inapplicable routes bypass RE auxiliaries", {
  local_mocked_bindings(panel_select = function(...) stop("Auxiliary should not run"),
    .package = "PanelSelect")
  fixture = probit_initial_fixture()
  control = fit_control(init = "re")
  complete = initial_parameters(fixture$prepared, fixture$info, NULL, numeric(),
    fit_control(init = "pooled"))
  local_test_seed(42); saved = .Random.seed
  expect_identical(initial_parameters(fixture$prepared, fixture$info, complete, numeric(), control), complete)
  expect_identical(.Random.seed, saved)
  for (fixture in list(probit_initial_fixture(cross = TRUE), probit_initial_fixture(selected = FALSE))) {
    expected = initial_parameters(fixture$prepared, fixture$info, NULL, numeric(),
      fit_control(init = "pooled"))
    expect_identical(initial_parameters(fixture$prepared, fixture$info, NULL, numeric(), control), expected)
  }
})

test_that("a supplied complete equation block is skipped and failure preserves pooled starts", {
  fixture = probit_initial_fixture()
  prepared = fixture$prepared; info = fixture$info
  pooled = initial_parameters(prepared, info, NULL, numeric(), fit_control(init = "pooled"))
  block = c(paste0("selection.", colnames(prepared$w)), "sd_re_selection")
  supplied = pooled[block]
  calls = 0L
  local_mocked_bindings(panel_select = function(...) {
    calls <<- calls + 1L
    stop("deliberate auxiliary failure")
  }, .package = "PanelSelect")
  expect_warning({ actual = initial_parameters(prepared, info, supplied, numeric(),
    fit_control(init = "re")) }, "outcome RE initialization failed")
  expect_identical(calls, 1L)
  expect_identical(actual, pooled)
})

test_that("auxiliary convergence matters but unused covariance and Q warnings do not", {
  fixture = probit_initial_fixture()
  pooled = initial_parameters(fixture$prepared, fixture$info, NULL, numeric(),
    fit_control(init = "pooled"))
  outcome = pooled[c(paste0("outcome.", colnames(fixture$prepared$x)), "sd_re_outcome")]
  converged = FALSE
  local_mocked_bindings(panel_select = function(outcome, selection, data, id, family, ...) {
    info = parameter_layout(prepare_data(outcome, NULL, data, id, family))
    coefficients = setNames(rep(.6, nrow(info)), info$name)
    if (!converged) warning("Optimizer did not report convergence: fixture", call. = FALSE)
    warning("The requested covariance estimate is unavailable; inspect the fit diagnostics", call. = FALSE)
    warning("Optimizer stopped with Q_OPG = 1.5", call. = FALSE)
    list(coefficients = coefficients, optimizer = list(converged = converged))
  }, .package = "PanelSelect")
  expect_warning({ fallback = initial_parameters(fixture$prepared, fixture$info, outcome,
    numeric(), fit_control(init = "re")) }, "selection RE initialization failed")
  expect_identical(fallback, pooled)
  converged = TRUE
  expect_silent({ improved = initial_parameters(fixture$prepared, fixture$info, outcome,
    numeric(), fit_control(init = "re")) })
  selection = startsWith(names(improved), "selection.") | names(improved) == "sd_re_selection"
  expect_equal(unname(improved[selection]), rep(.6, sum(selection)))
  expect_identical(improved[!selection], pooled[!selection])
})

test_that("native ProbitRE starts agree with separate small-panel auxiliary fits", {
  fixture = probit_initial_fixture()
  data = fixture$data
  prepared = prepare_data(y ~ x + offset(off), d ~ w + offset(off), data, "id", "probit")
  info = parameter_layout(prepared)
  for (method in c("gh", "agh")) {
    control = fit_control(panel_control(quadrature = method, nodes = c(3, 5)), init = "re")
    initial = initial_parameters(prepared, info, NULL, numeric(), control)
    pooled = initial_parameters(prepared, info, NULL, numeric(), fit_control(init = "pooled"))
    selection_names = c("selection.(Intercept)", "selection.w", "sd_re_selection")
    outcome_names = c("outcome.(Intercept)", "outcome.x", "sd_re_outcome")
    selection_start = setNames(pooled[selection_names], c("outcome.(Intercept)", "outcome.w", "sd_re_outcome"))
    selection = panel_select(d ~ w + offset(off), NULL, data, "id", "probit",
      start = selection_start, control = panel_control(quadrature = method, nodes = 3), vcov = "opg")
    # The small selected-outcome sample stops with a sizeable Q_OPG.
    outcome = with_warnings(panel_select(y ~ x + offset(off), NULL, data[data$d == 1, ], "id", "probit",
      start = pooled[outcome_names], control = panel_control(quadrature = method, nodes = 5), vcov = "opg"))
    expect_length(outcome$warnings, 1L)
    expect_match(outcome$warnings, 'Q_OPG =')
    outcome = outcome$value
    expect_true(selection$optimizer$converged)
    expect_true(outcome$optimizer$converged)
    expect_equal(unname(initial[selection_names]), unname(coef(selection)), tolerance = 1e-10)
    expect_equal(initial[outcome_names], coef(outcome), tolerance = 1e-10)
    expect_true(all(is.finite(initial)))
  }
})

test_that("offset-only probit auxiliaries preserve zero fixed and supplied scales", {
  fixture = probit_initial_fixture(no_intercept = TRUE)
  info = fixture$info
  fixed = resolve_fixed(c(sd_re_selection = 0), info, fixture$prepared)
  info$fixed = info$name %in% names(fixed)
  start = c(sd_re_outcome = .45)
  actual = initial_parameters(fixture$prepared, info, start, fixed, fit_control(init = "re"))
  expect_identical(actual[names(fixed)], fixed)
  expect_identical(actual[names(start)], start)
  expect_false(any(startsWith(names(actual), "selection.") | startsWith(names(actual), "outcome.")))
  # Releasing the outcome scale actually fits a one-parameter offset-only auxiliary.
  actual = initial_parameters(fixture$prepared, info, NULL, fixed,
    fit_control(panel_control(nodes = 3), init = "re"))
  expect_true(is.finite(actual[["sd_re_outcome"]]) && actual[["sd_re_outcome"]] > 0)
  expect_identical(actual[names(fixed)], fixed)
})
