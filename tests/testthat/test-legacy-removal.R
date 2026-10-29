# Removed interfaces include the PanelSelect 1.x and PanelCount wrappers,
# their predictions, check_binary_quadrature(), the checked binary alias, frozen
# adaptive geometry and legacy optimizer options. Calls name the replacement.
# Historical objects in fixtures/saved-legacy-settings.R remain readable:
# default-valued compatibility fields are dropped; fits with a removed setting
# display but stop where they would use it, with instructions for refitting.

legacy = dget(test_path('fixtures', 'saved-legacy-settings.R'))

test_that('only the modern functions are exported', {
  expect_identical(sort(getNamespaceExports('PanelSelect')), sort(c('add_information',
    'check_quadrature', 'check_stationarity', 'panel_control', 'panel_select', 'random_effects')))
  removed = c('probitRE_linearRE', 'probitRE_probitRE', 'probitRE_PoissonRE', 'probitRE_PLNRE',
    'ProbitRE_PoissonRE', 'ProbitRE_PLNRE', 'ProbitRE', 'PoissonRE', 'PLN_RE',
    'predict_probitRE_linearRE', 'predict_probitRE_probitRE', 'predict_ProbitRE_PoissonRE',
    'predict_ProbitRE_PLNRE', 'check_binary_quadrature')
  for (name in removed) expect_false(exists(name, asNamespace('PanelSelect'), inherits = FALSE), label = name)
})

test_that('the removed compatibility settings stop with their replacements', {
  expect_error(panel_control(binary_integration = 'checked'),
    'binary_integration was removed; use quadrature = "gl"', fixed = TRUE)
  expect_error(panel_control(checked = list(orders = c(10, 20))),
    'checked was removed with binary_integration; use quadrature = "gl" and quadrature_options',
    fixed = TRUE)
  expect_error(panel_control(fixed_geometry = list()),
    'fixed_geometry was removed; use quadrature = "agh"', fixed = TRUE)
  expect_false(any(c('fixed_geometry', 'binary_integration', 'checked') %in% names(panel_control())))
})

test_that('controls saved with default-valued compatibility fields read as current controls', {
  expect_identical(current_control(legacy$controls$default), panel_control(quadrature = 'agh'))
  expect_identical(current_control(legacy$controls$gl), panel_control(quadrature = 'gl'))
  expect_identical(validate_panel_control(legacy$controls$gl), panel_control(quadrature = 'gl'))
  old = legacy$fits$modern
  fit = panel_select(count ~ x, NULL, legacy$data$count, 'id', 'poisson',
    control = panel_control(quadrature = 'gh', nodes = 6L))
  expected_control = fit$control; expected_control$init = old$control$init
  expect_identical(current_control(old$control), expected_control)
  # A historical fit and a fresh optimization can differ by floating-point roundoff.
  expect_equal(coef(old), coef(fit), tolerance = 1e-8)
  expect_equal(predict(old, legacy$data$count), predict(fit, legacy$data$count), tolerance = 1e-8)
  expect_equal(coef(update(old, init = 're')), coef(update(fit, init = 're')), tolerance = 1e-8)
})

test_that('fits saved with the checked alias or frozen geometry stop where they would integrate', {
  checked = legacy$fits$checked
  hint = 'refit with update(fit, control = panel_control(quadrature = "gl"))'
  expect_no_error(utils::capture.output(print(summary(checked))))
  expect_true(all(is.finite(vcov(checked))))
  expect_error(update(checked), paste('binary_integration = "checked" was removed;', hint),
    fixed = TRUE)
  expect_error(check_quadrature(checked), hint, fixed = TRUE)
  expect_error(add_information(checked), hint, fixed = TRUE)
  # The checked alias evaluated GL with orders 20 and 10, the GL defaults: the
  # suggested refit has the same objective and optimum.
  same = panel_select(binary ~ x, d ~ x + z, legacy$data$selection, 'id', 'probit', fixed = coef(checked),
    control = panel_control(quadrature = 'gl'))
  expect_equal(same$loglik, checked$loglik, tolerance = 1e-12)
  refit = update(checked, control = panel_control(quadrature = 'gl'), init = 'pooled')
  expect_identical(refit$optimizer$status, 'converged')
  expect_lt(abs(refit$loglik - checked$loglik), 1e-6)

  frozen = legacy$fits$frozen
  hint = 'refit with update(fit, control = panel_control(quadrature = "agh"))'
  expect_no_error(utils::capture.output(print(summary(frozen))))
  expect_error(update(frozen), paste('frozen adaptive geometry (fixed_geometry) was removed;', hint),
    fixed = TRUE)
  expect_error(check_quadrature(frozen), hint, fixed = TRUE)
  expect_error(random_effects(frozen), hint, fixed = TRUE)
  expect_error(predict(frozen, legacy$data$count, conditioning = 'history', history = legacy$data$count),
    hint, fixed = TRUE)
})

test_that('a fit saved with legacy optimizer options evaluates but stops on refits', {
  fit = legacy$fits$legacy_options
  message = paste('legacy optimizer options (REPORT) were removed; refit with',
    'update(fit, control = panel_control()), adding the settings to keep')
  expect_no_error(utils::capture.output(print(summary(fit))))
  expect_s3_class(check_quadrature(fit), 'panelselect_quadrature_check')
  expect_error(update(fit), message, fixed = TRUE)
  expect_error(check_quadrature(fit, refit = TRUE), message, fixed = TRUE)
  refit = update(fit, control = panel_control(quadrature = 'gh', nodes = 4L, maxit = 0L), init = 'pooled')
  expect_null(attr(refit$control, 'legacy_optimizer_control'))
  expect_identical(refit$control$quadrature, 'gh')
})
