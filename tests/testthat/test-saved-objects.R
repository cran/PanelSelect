# Actual historical objects in fixtures/saved-vocabulary.R are read in the
# current vocabulary. Saved controls normalize to current controls. Each fit is
# rebuilt from its stored data with an equivalent current call; the historical
# components then replace the corresponding current components. Normalization
# must preserve stored results exactly; independent refits need numerical agreement,
# not identical optimization paths. Same-fit vocabulary tests check exact method parity.

saved = dget(test_path('fixtures', 'saved-vocabulary.R'))
saved_budget = list(seconds_per_panel = 1e4, seconds_total = 1e4)

# Timings differ by run and are excluded from numerical comparisons.
timeless = function(x) {
  if (!is.list(x)) return(x)
  x[names(x) %in% c('time', 'elapsed_seconds', 'total_seconds', 'panel_seconds', 'active_row',
    'inner_time')] = NULL
  for (i in seq_along(x)) if (is.list(x[[i]])) x[[i]] = timeless(x[[i]])
  x
}

# Ignore translated labels, retaining every numeric value, type and array shape.
numeric_payload = function(x) {
  if (is.numeric(x)) return(unname(x))
  if (is.list(x)) return(unname(lapply(x, numeric_payload)))
  NULL
}

# Explicit AGH preserves the default used when the historical fixtures were saved.
rebuild = function(name) {
  data = saved$data
  suppressWarnings(switch(name,
    linear_gh = panel_select(y ~ x, d ~ w, data$selection, 'id', 'gaussian',
      control = panel_control(quadrature = 'agh', nodes = 4L, maxit = 200L), init = 'pooled'),
    count_gl = panel_select(count ~ x, NULL, data$count, 'id', 'poisson',
      control = panel_control(quadrature = 'gl', nodes = 12L, quadrature_options = saved_budget,
        maxit = 200L)),
    pln_gl_observation = panel_select(count ~ x, NULL, data$count, 'id', 'pln',
      control = panel_control(quadrature = 'gh', nodes = 6L, observation_quadrature = 'gl',
        observation_nodes = 20L, observation_options = c(saved_budget, list(max_cells = 16384L)),
        maxit = 200L)),
    cross_linear = panel_select(y ~ x, d ~ w, data$cross, NULL, 'gaussian',
      control = panel_control(quadrature = 'agh', maxit = 200L)),
    poisson_fixed = panel_select(count ~ x, NULL, data$poisson, 'id', 'poisson',
      fixed = c(sd_re_outcome = 0), control = panel_control(quadrature = 'agh', maxit = 200L)),
    # The fit of the removed 1.x call probitRE_linearRE(d ~ w, y ~ x, 'id',
    # data$selection, H = 4, use.optim = TRUE, verbose = -1).
    legacy_linear = panel_select(y ~ x, d ~ w, data$selection, 'id', 'gaussian',
      start = c(cor_re = 0, cor_error = 0, sd_re_outcome = 1, sd_error = 1),
      control = panel_control(quadrature = 'gh', nodes = 4L, maxit = 100L), init = 'pooled')))
}

# The rebuilt fit with the historical components in their places.
as_saved = function(fit, parts) {
  old = unclass(fit)
  names(old)[names(old) == 'loglik_by_unit'] = 'loglik_by_panel'
  for (name in setdiff(names(parts), 'outcome_scale')) old[name] = list(parts[[name]])
  if (!is.null(parts$outcome_scale)) old$prepared$outcome_scale = parts$outcome_scale
  structure(old, class = class(fit))
}

test_that('historical controls normalize to the controls the current interface builds', {
  current = list(
    default = fit_control(panel_control(quadrature = 'agh', maxit = 200L)),
    fixed_rule = fit_control(panel_control(quadrature = 'gh', nodes = 12L, maxit = 200L)),
    split_options = fit_control(panel_control(quadrature = 'gl', nodes = 16L, maxit = 200L,
      quadrature_options = list(panel_loglik_tol = 1e-8, max_cells = 2048L))),
    inner_adaptive = fit_control(panel_control(quadrature = 'agh', observation_quadrature = 'agh', observation_nodes = 12L,
      maxit = 200L)),
    inner_split = fit_control(panel_control(quadrature = 'agh', observation_quadrature = 'gl', observation_nodes = 15L,
      observation_options = list(panel_score_tol = 1e-6), maxit = 200L)),
    no_reduction = fit_control(panel_control(quadrature = 'agh', exact_reductions = FALSE, maxit = 200L)),
    no_stationarity = fit_control(panel_control(quadrature = 'agh', stationarity_warning = NULL, maxit = 200L)),
    parscale = fit_control(panel_control(quadrature = 'agh', parscale = c(a = 1, b = 2), maxit = 200L)),
    settings = fit_control(panel_control(quadrature = 'agh', maxit = 200L), na_action = 'omit', init = 'pooled',
      vcov = 'observed'),
    maxlik = fit_control(panel_control(quadrature = 'agh', optimizer = 'maxLik', method = 'BHHH', maxit = 50L,
      reltol = 1e-10)))
  expect_identical(setdiff(names(saved$controls), c('checked', 'gk_options')), names(current))
  for (name in names(current)) {
    # Display preserves the stored initializer, including ProbitRE.
    current[[name]]$init = saved$controls[[name]]$init
    expect_identical(current_control(saved$controls[[name]]), current[[name]], label = name)
  }
  # A saved GK control is renamed for display and
  # stops where it would be revalidated, naming the GL refit.
  gk = current_control(saved$controls$gk_options)
  expect_identical(gk$quadrature, 'gk')
  expect_identical(gk$quadrature_options[c('total_score_tol', 'max_bound')],
    list(total_score_tol = 1e-3, max_bound = 32))
  expect_error(validate_panel_control(saved$controls$gk_options),
    'GK integration was removed; refit with update(fit, control = panel_control(quadrature = "gl"))',
    fixed = TRUE)
  # The removed checked binary alias receives the same treatment.
  checked = current_control(saved$controls$checked)
  expect_identical(checked$binary_integration, 'checked')
  expect_identical(checked$checked$orders, c(10L, 20L))
  expect_error(validate_panel_control(saved$controls$checked),
    'binary_integration = "checked" was removed; refit with update(fit, control = panel_control(quadrature = "gl"))',
    fixed = TRUE)
})

test_that('historical fits preserve stored results and agree numerically with current fits', {
  for (name in names(saved$fits)) {
    fit = rebuild(name)
    old = as_saved(fit, saved$fits[[name]])
    normalized = current_fit(old)
    normalized$call = fit$call
    expect_identical(normalized$control$init, saved$fits[[name]]$control$init, label = name)
    fit$control$init = saved$fits[[name]]$control$init
    if (name == 'legacy_linear') {
      # The name map of the removed wrapper is kept, in current parameter names.
      expect_identical(normalized$legacy, list(style = 'gaussian', parameter_names = c(
        'selection.(Intercept)' = 'probit.(Intercept)', selection.w = 'probit.w',
        'outcome.(Intercept)' = 'linear.(Intercept)', outcome.x = 'linear.x', sd_re_selection = 'delta',
        sd_re_outcome = 'lambda', sd_error = 'sigma', cor_re = 'rho', cor_error = 'tau')))
      normalized$legacy = NULL
    }
    # Historical effort counts and diagnostic probes belong to their original run.
    # Only stationarity changes convention, from score_norm/quadratic to Q.
    fields = c('coefficients', 'covariance', 'loglik', 'score', 'optimizer',
      'integration', 'fixed', 'parameter_info')
    stored = old[fields]; current = normalized[fields]
    stored$optimizer$stationarity = current$optimizer$stationarity = NULL
    expect_identical(numeric_payload(current), numeric_payload(stored), label = name)
    expect_identical(normalized$loglik_by_unit, old$loglik_by_panel, label = name)
    old_mesh = old$control$integration_mesh; new_mesh = normalized$control$integration_mesh
    if (!is.null(old_mesh)) {
      old_mesh$binding = new_mesh$binding = NULL
      expect_identical(numeric_payload(new_mesh), numeric_payload(old_mesh), label = name)
    }
    gradient = normalized$optimizer$gradient
    quadratic = sum(gradient * (normalized$optimizer$covariance$opg %*% gradient))
    expect_equal(normalized$optimizer$stationarity$Q_OPG, quadratic, tolerance = 1e-12,
      label = name)
    expect_identical(normalized$optimizer$stationarity$tolerance,
      normalized$control$stationarity_warning, label = name)

    # Structure and settings remain exact; floating-point refits use a small tolerance.
    expect_identical(class(normalized), class(fit), label = name)
    expect_identical(names(normalized), names(fit), label = name)
    fields = c('model_id', 'structure', 'parameter_info', 'fixed', 'prepared')
    expect_identical(normalized[fields], fit[fields], label = name)
    control_fields = setdiff(names(fit$control), 'integration_mesh')
    expect_identical(normalized$control[control_fields], fit$control[control_fields], label = name)
    fields = c('coefficients', 'covariance', 'loglik', 'loglik_by_unit', 'score')
    expect_equal(normalized[fields], fit[fields], tolerance = 1e-8, label = name)
  }
})

test_that('every method treats a historical fit as the current fit', {
  rows = saved$data$selection[1:4, ]
  for (name in c('linear_gh', 'count_gl', 'legacy_linear')) {
    fit = rebuild(name)
    old = as_saved(fit, saved$fits[[name]])
    old$call = fit$call
    targets = if (name == 'count_gl') saved$data$count[1:3, ] else rows
    expect_equal(coef(old), coef(fit), tolerance = 1e-8, label = name)
    expect_equal(vcov(old), vcov(fit), tolerance = 1e-8, label = name)
    expect_equal(vcov(old, 'sandwich'), vcov(fit, 'sandwich'), tolerance = 1e-8, label = name)
    expect_equal(confint(old), confint(fit), tolerance = 1e-8, label = name)
    expect_equal(predict(old, targets, se.fit = TRUE), predict(fit, targets, se.fit = TRUE),
      tolerance = 1e-8, label = name)
    expect_equal(residuals(old), residuals(fit), tolerance = 1e-8, label = name)
    expect_equal(summary(old)$coefficients, summary(fit)$coefficients, tolerance = 1e-8,
      label = name)
    # Printed effort counts describe the saved run, not a new optimizer trajectory.
    counts = old$optimizer$counts
    printed = paste(utils::capture.output(print(summary(old))), collapse = '\n')
    expect_match(printed, paste(counts[[1L]], 'function and', counts[[2L]],
      'gradient evaluations'), fixed = TRUE)
    expect_equal(timeless(check_quadrature(old)), timeless(check_quadrature(fit)),
      tolerance = 1e-8, label = name)
    # Refitting a stored ProbitRE policy requires an explicit current initializer.
    old_call = update(old, init = fit$control$init, evaluate = FALSE)
    new_call = update(fit, init = fit$control$init, evaluate = FALSE)
    expect_identical(old_call$start, coef(old), label = name)
    expect_equal(old_call$start, new_call$start, tolerance = 1e-8, label = name)
    old_call$start = new_call$start
    expect_identical(old_call, new_call, label = name)
    expect_equal(coef(update(old, init = fit$control$init)),
      coef(update(fit, init = fit$control$init)), tolerance = 1e-8, label = name)
  }
  # A retained staged mesh binds after current_control() alone, as validation uses it.
  fit = rebuild('count_gl'); old = as_saved(fit, saved$fits$count_gl)
  expect_identical(timeless(evaluate_model(coef(fit), fit$prepared, current_control(old$control), fit$fixed)),
    timeless(evaluate_model(coef(fit), fit$prepared, fit$control, fit$fixed)))
})

test_that('a historical summary and quadrature check still print', {
  expect_no_error(utils::capture.output(print(saved$summary)))
  fit = rebuild('linear_gh')
  expect_match(utils::capture.output(print(saved$check))[1L], '^Quadrature check [(]higher_order[)] of: ')
  # The compared rules and changes match; the historical check also filled the observation
  # columns of models without a PLN observation integral.
  old = current_check(saved$check); new = check_quadrature(fit)
  columns = c('rule', 'nodes', 'loglik_change', 'max_panel_change', 'max_score_change',
    'max_panel_score_change', 'within_tolerance')
  exact = c('rule', 'nodes', 'within_tolerance')
  expect_identical(old$comparison[exact], new$comparison[exact])
  expect_equal(old$comparison[setdiff(columns, exact)], new$comparison[setdiff(columns, exact)],
    tolerance = 1e-8)
  expect_true(all(is.na(new$comparison$observation_nodes)))
})
