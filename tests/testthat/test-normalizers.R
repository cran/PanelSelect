# Saved controls and fits are read in the current vocabulary by
# current_control() and current_fit(); controls and fits made now pass through
# unchanged, and a control that mixes the two vocabularies stops.

normalizer_data = function(seed = 1501, panels = 24L, waves = 3L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  e = stats::rnorm(n); epsilon = .3*e + sqrt(1 - .3^2)*stats::rnorm(n)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + e > 0)
  data$y = ifelse(data$d == 1, .5 + .4*data$x + .6*v[data$id] + .5*epsilon, NA)
  data$count = stats::rpois(n, exp(.1 + .3*data$x + .8*u[data$id]))
  data
}

test_that('current controls and fits pass through the normalizers unchanged', {
  data = normalizer_data()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'linear',
    control = panel_control(quadrature = 'agh', nodes = 4L),
    init = 'pooled')
  expect_identical(current_fit(fit), fit)
  expect_identical(current_control(fit$control), fit$control)
  control = panel_control(quadrature = 'gl', quadrature_options = list(panel_loglik_tol = 1e-8))
  expect_identical(current_control(control), control)
  once = validate_panel_control(control)
  expect_identical(validate_panel_control(once), once)
  # Normalizing an earlier fit twice changes nothing more.
  old = old_vocabulary(fit)
  expect_identical(current_fit(current_fit(old)), current_fit(old))
})

test_that('earlier names added to a current control stop instead of being overwritten', {
  control = panel_control(quadrature = 'gl', quadrature_options = list(bound = 10))
  control$inner_nodes = 40
  expect_error(validate_panel_control(control),
    'control mixes earlier and current names: inner_nodes is now observation_nodes', fixed = TRUE)
  stored = panel_control()
  stored$reduction = FALSE
  expect_error(validate_panel_control(stored), 'reduction is now exact_reductions', fixed = TRUE)
  stored = panel_control()
  stored$split = list(loglik_tol = 1e-8)
  expect_error(validate_panel_control(stored), 'split is now quadrature_options', fixed = TRUE)
  rule = panel_control()
  rule$quadrature = 'adaptive'
  expect_error(validate_panel_control(rule), 'quadrature = "adaptive" is now "agh"', fixed = TRUE)
  data = normalizer_data()
  expect_error(panel_select(count ~ x, NULL, data, 'id', 'poisson', control = rule),
    'control mixes earlier and current names')
})

test_that('labels derived from the data keep their names in saved fits', {
  data = normalizer_data()
  # Panel IDs spelled like the earlier parameter names.
  ids = paste0('p', seq_len(24L)); ids[1:2] = c('sd_outcome', 'sd_selection')
  data$id = ids[data$id]
  fit = panel_select(count ~ x, NULL, data, 'id', 'poisson', control = panel_control(nodes = 4L))
  normalized = current_fit(old_vocabulary(fit))
  expect_identical(names(normalized$loglik_by_unit), names(fit$loglik_by_unit))
  expect_identical(normalized$data, fit$data)
  expect_identical(rownames(normalized$score), rownames(fit$score))
  expect_identical(normalized$coefficients, fit$coefficients)
  # A regressor named like a parameter: the removed count wrappers
  # stored a name map and estimates labelled by column name.
  data$sd_outcome = data$x
  fit = panel_select(count ~ sd_outcome, NULL, data, 'id', 'poisson', control = panel_control(nodes = 4L))
  fit$legacy = list(style = 'poisson', parameter_names = c('outcome.(Intercept)' = '(Intercept)',
    outcome.sd_outcome = 'sd_outcome', sd_re_outcome = 'sigma'))
  fit$par = setNames(unname(fit$coefficients), unname(fit$legacy$parameter_names))
  old = old_vocabulary(fit)
  expect_true('sd_outcome' %in% names(old$legacy$parameter_names))
  normalized = current_fit(old)
  expect_identical(names(normalized$par), names(fit$par))
  expect_identical(normalized$legacy, fit$legacy)
  expect_identical(normalized$coefficients, fit$coefficients)
})

test_that('a saved control keeps its attributes and prints in the current names', {
  data = normalizer_data()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'linear',
    control = panel_control(quadrature = 'agh', nodes = 4L, reltol = 1e-9), init = 'pooled')
  old = old_vocabulary(fit)$control
  expect_true(attr(old, 'reltol_supplied'))
  normalized = current_control(old)
  expect_identical(attributes(normalized), attributes(fit$control))
  expect_identical(utils::capture.output(print(old)), utils::capture.output(print(fit$control)))
  # Empty option lists come back as list(), as panel_control() stores them.
  expect_identical(normalized$quadrature_options, list())
})

test_that('saved moving-cut GL diagnostics report their options in the current names', {
  data = data.frame(id = rep(1:3, each = 2), d = c(1, 1, 1, 0, 0, 0), y = c(1, 1, 0, NA, NA, NA))
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0,
    sd_re_selection = 1, sd_re_outcome = 1, cor_re = 0, cor_error = 0)
  fit = panel_select(y ~ 1, d ~ 1, data, 'id', 'probit', fixed = theta,
    control = panel_control(quadrature = 'gl', nodes = 20))
  current = c('panel_loglik_tol', 'panel_score_tol', 'total_loglik_tol', 'total_score_tol', 'max_cells')
  expect_true(all(current %in% names(fit$integration$options)))
  # Historical controls reported the options under their earlier names.
  old = old_vocabulary(fit)
  earlier = c('loglik_tol', 'score_tol', 'aggregate_loglik_tol', 'aggregate_score_tol', 'max_intervals')
  at = match(current, names(old$integration$options))
  names(old$integration$options)[at] = earlier
  expect_identical(current_fit(old)$integration, fit$integration)
  expect_identical(summary(old)$integration$options, summary(fit)$integration$options)
})
