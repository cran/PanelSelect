# Public names: integration rules agh, gh and gl; one options
# list per integral; panel_*/total_* tolerances and max_cells; exact_reductions
# and stationarity_warning; check_quadrature(compare =); prediction targets
# selection_probability and selected_mean; person-effect SDs sd_re_selection
# and sd_re_outcome, their correlation cor_re; and the fit component loglik_by_unit.
# Earlier names in new calls stop with their replacements.

# Person effects correlated 0.4 and errors correlated 0.3: interior values.
names_data = function(seed = 1301, panels = 30L, waves = 4L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  e = stats::rnorm(n); epsilon = .3*e + sqrt(1 - .3^2)*stats::rnorm(n)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + e > 0)
  data$y = ifelse(data$d == 1, .5 + .4*data$x + .6*v[data$id] + .5*epsilon, NA)
  data$count = stats::rpois(n, exp(.1 + .3*data$x + .5*u[data$id]))
  data
}

names_fit = function(data = names_data(), ...)
  panel_select(y ~ x, d ~ w, data, 'id', 'linear', control = panel_control(nodes = 4L, ...),
    init = 'pooled')

test_that('earlier rule names stop with their replacements', {
  replacements = c(adaptive = 'agh', fixed = 'gh', split = 'gl')
  for (old in names(replacements)) {
    expect_error(panel_control(quadrature = old),
      sprintf('quadrature = "%s" is now "%s"', old, replacements[[old]]), fixed = TRUE)
    expect_error(panel_control(observation_quadrature = old),
      sprintf('observation_quadrature = "%s" is now "%s"', old, replacements[[old]]), fixed = TRUE)
  }
})

test_that('earlier argument and option-list names stop with their replacements', {
  renamed = c(inner_quadrature = 'observation_quadrature', inner_nodes = 'observation_nodes',
    reduction = 'exact_reductions', stationarity_tol = 'stationarity_warning',
    split = 'quadrature_options', gk = 'quadrature_options',
    inner_split = 'observation_options', inner_gk = 'observation_options')
  for (old in names(renamed))
    expect_error(do.call(panel_control, stats::setNames(list(NULL), old)),
      paste(old, 'is now', renamed[[old]]), fixed = TRUE)
  options = c(loglik_tol = 'panel_loglik_tol', score_tol = 'panel_score_tol',
    aggregate_loglik_tol = 'total_loglik_tol', aggregate_score_tol = 'total_score_tol',
    max_intervals = 'max_cells')
  for (old in names(options))
    expect_error(panel_control(quadrature = 'gl',
      quadrature_options = stats::setNames(list(1e-3), old)),
      paste('quadrature_options entry', old, 'is now', options[[old]]), fixed = TRUE)
  expect_error(panel_control(observation_quadrature = 'gl',
    observation_options = list(max_intervals = 10)),
    'observation_options entry max_intervals is now max_cells', fixed = TRUE)
})

test_that('the new defaults and the options of each rule are explicit', {
  control = panel_control()
  expect_identical(control$quadrature, 'gh')
  expect_identical(control$observation_quadrature, 'gh')
  expect_identical(control$quadrature_options, list())
  expect_identical(control$observation_options, list())
  expect_true(control$exact_reductions)
  expect_identical(control$stationarity_warning, .1)
  expect_identical(control$maxit, 500L)
  expect_false(any(c('split', 'gk', 'inner_split', 'inner_gk', 'inner_quadrature', 'inner_nodes',
    'reduction', 'stationarity_tol', 'na_action', 'init', 'vcov') %in% names(control)))
  gl = list(reference_nodes = 10L, bound = 12, panel_loglik_tol = 1e-7, panel_score_tol = 1e-5,
    tail_tol = 1e-9, total_loglik_tol = 1e-5, total_score_tol = 1e-4, information_tol = 2e-3,
    seconds_per_panel = 60, seconds_total = 300, max_nodes = 100000000L, max_cells = 4096L,
    batch_nodes = 8192L)
  expect_identical(panel_control(quadrature = 'gl')$quadrature_options, gl)
  observation = gl; observation$batch_nodes = 128L
  expect_identical(panel_control(observation_quadrature = 'gl')$observation_options, observation)
  # Options are checked against the chosen rule.
  expect_error(panel_control(quadrature = 'gh', quadrature_options = list(bound = 8)),
    "quadrature_options require quadrature = 'gl'", fixed = TRUE)
  expect_error(panel_control(observation_options = list(bound = 8)),
    "observation_options require observation_quadrature = 'gl'", fixed = TRUE)
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(bounds = 8)),
    'Unknown quadrature_options entry for the gl rule: bounds', fixed = TRUE)
  expect_error(panel_control(exact_reductions = NA), 'exact_reductions must be TRUE or FALSE')
  expect_error(panel_control(stationarity_warning = -1), 'stationarity_warning must be NULL')
})

test_that('rule names from another vocabulary cannot fall through to fixed Gauss-Hermite', {
  expect_error(prepare_rule_control(list(quadrature = 'adaptive', nodes = 5L,
    observation_quadrature = 'gh', observation_nodes = 5L)), 'unknown integration rule')
  expect_error(prepare_rule_control(list(quadrature = 'gh', nodes = 5L,
    observation_quadrature = 'fixed', observation_nodes = 5L)), 'unknown integration rule')
  expect_error(staged_options(list(quadrature = 'gl', observation_quadrature = 'split'), inner = TRUE),
    'no staged observation rule named split')
  # Stored controls in the earlier vocabulary are renamed before evaluation.
  fixture = core_fixture('selection_pln')
  current = core_control('agh', nodes = 3L, inner = 5L)
  earlier = list(quadrature = 'adaptive', nodes = 3L, inner_nodes = 5L, reduction = TRUE,
    mode_maxit = 60L, mode_tol = 1e-10, fixed_geometry = NULL)
  expect_identical(evaluate_model(fixture$theta, fixture$prepared, earlier),
    evaluate_model(fixture$theta, fixture$prepared, current))
})

test_that('check_quadrature() compares with compare = and names observation orders', {
  fit = names_fit()
  check = check_quadrature(fit)
  expect_identical(check$compare, 'higher_order')
  expect_null(check$method)
  expect_true(all(c('observation_nodes', 'observation_quadrature') %in% names(check$comparison)))
  # A model without a PLN observation integral has no observation columns to fill.
  expect_true(all(is.na(check$comparison$observation_nodes)))
  expect_true(all(is.na(check$comparison$observation_quadrature)))
  # The comparison is the change from the fitted rule at the estimates.
  refined = fit$control; refined$nodes = 9L
  expected = sum(evaluate_model(coef(fit), fit$prepared, refined, fit$fixed)$loglik -
    evaluate_model(coef(fit), fit$prepared, fit$control, fit$fixed)$loglik)
  expect_identical(check$comparison$loglik_change[check$comparison$rule == 'person_effect'], expected)
  expect_output(print(check), 'Quadrature check \\(higher_order\\) of: ')
  expect_error(check_quadrature(fit, method = 'refine'),
    'method is now compare; method = "refine" is now compare = "higher_order"', fixed = TRUE)
  expect_error(check_quadrature(fit, inner_nodes = 30), 'inner_nodes is now observation_nodes')
  expect_error(check_quadrature(fit, gk = list()), 'gk is now quadrature_options')
  expect_error(check_quadrature(fit, compare = 'refine'), 'compare = "refine" is now "higher_order"',
    fixed = TRUE)
  expect_error(check_quadrature(fit, other = 1), 'does not use the argument')
  # Requested options are checked before integration: an earlier or unknown
  # name stops rather than reporting an integration failure.
  expect_error(check_quadrature(fit, compare = 'gl', quadrature_options = list(max_intervals = 10)),
    'quadrature_options entry max_intervals is now max_cells', fixed = TRUE)
  expect_error(check_quadrature(fit, compare = 'gl', quadrature_options = list(unknown = 10)),
    'Unknown quadrature_options entry for the gl rule: unknown', fixed = TRUE)
  staged = panel_select(count ~ x, NULL, names_data(panels = 6L), 'id', 'poisson',
    control = staged_test_budget(panel_control(quadrature = 'gl', nodes = 12L)))
  expect_error(check_quadrature(staged, nodes = 300), 'GL nodes must be one integer from 3 through 256',
    fixed = TRUE)
})

test_that('a fit with no person effect left reports no person-effect integral', {
  cross = names_data(); cross$id = NULL
  fit = panel_select(y ~ x, d ~ w, cross, NULL, 'linear', control = panel_control(quadrature = 'gl'))
  expect_identical(summary(fit)$integration_description, 'None needed (closed-form likelihood)')
})

test_that('prediction targets are selection_probability and selected_mean', {
  data = names_data()
  fit = names_fit(data)
  rows = data[1:4, ]
  probability = predict(fit, rows, type = 'selection_probability')
  joint = predict(fit, rows, type = 'joint')
  expect_equal(predict(fit, rows, type = 'selected_mean'), joint/probability, tolerance = 1e-12)
  expect_true(all(probability > 0 & probability < 1))
  expect_error(predict(fit, rows, type = 'selection'),
    'type = "selection" is now "selection_probability"', fixed = TRUE)
  expect_error(predict(fit, rows, type = 'selected'),
    'type = "selected" is now "selected_mean"', fixed = TRUE)
  # residuals() subtracts the selected mean by default.
  expect_identical(residuals(fit), residuals(fit, type = 'selected_mean'))
  expect_identical(predict(fit, rows, type = 'link', equation = 'selection'),
    drop(cbind(1, rows$w) %*% coef(fit)[c('selection.(Intercept)', 'selection.w')]) + 0,
    ignore_attr = TRUE)
})

test_that('person-effect parameters use re names and reject earlier names in new calls', {
  data = names_data()
  fit = names_fit(data)
  expect_identical(names(coef(fit)), c('selection.(Intercept)', 'selection.w',
    'outcome.(Intercept)', 'outcome.x', 'sd_re_selection', 'sd_re_outcome', 'sd_error',
    'cor_re', 'cor_error'))
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'linear', start = c(sd_selection = 1)),
    'sd_selection is now sd_re_selection (in start)', fixed = TRUE)
  expect_error(panel_select(count ~ x, NULL, data, 'id', 'poisson', fixed = c(sd_outcome = 0)),
    'sd_outcome is now sd_re_outcome (in fixed)', fixed = TRUE)
  renamed = c(sd_id_selection = 'sd_re_selection', sd_id_outcome = 'sd_re_outcome', cor_id = 'cor_re')
  for (old in names(renamed)) {
    value = stats::setNames(coef(fit)[[renamed[[old]]]], old)
    for (argument in c('start', 'fixed')) {
      supplied = stats::setNames(list(value), argument)
      message = paste0(old, ' is now ', renamed[[old]], ' (in ', argument, ')')
      expect_error(do.call(panel_select, c(list(y ~ x, d ~ w, data, 'id', 'linear'), supplied)),
        message, fixed = TRUE)
      expect_error(do.call(update, c(list(fit), supplied)), message, fixed = TRUE)
      expect_error(do.call(update, c(list(fit, evaluate = FALSE), supplied)), message, fixed = TRUE)
    }
  }
  natural = coef(fit)[unname(renamed)]
  updated_call = update(fit, start = natural, fixed = natural, evaluate = FALSE)
  expect_identical(updated_call$start[names(natural)], natural)
  expect_identical(updated_call$fixed, natural)
  effects = random_effects(fit)
  expect_identical(names(effects), c('id', 'selection_mean', 'selection_sd', 'outcome_mean',
    'outcome_sd'))
  expect_true(all(effects$selection_sd > 0 & effects$selection_sd < coef(fit)[['sd_re_selection']]))
})

test_that('loglik_by_unit is named by panel ID or by row name', {
  data = names_data()
  fit = names_fit(data)
  expect_null(fit$loglik_by_panel)
  expect_identical(names(fit$loglik_by_unit), as.character(unique(data$id)))
  expect_equal(sum(fit$loglik_by_unit), fit$loglik)
  cross = data; cross$id = NULL
  cross_fit = panel_select(y ~ x, d ~ w, cross, NULL, 'linear')
  expect_identical(names(cross_fit$loglik_by_unit), rownames(cross))
})

test_that('the summary legend lists the structural parameters of the model only', {
  data = names_data()
  legend = function(fit) names(summary(fit)$parameter_legend)
  expect_identical(legend(names_fit(data)),
    c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error'))
  expect_identical(legend(panel_select(count ~ x, NULL, data, 'id', 'poisson')), 'sd_re_outcome')
  cross = data; cross$id = NULL
  expect_identical(legend(panel_select(y ~ x, d ~ w, cross, NULL, 'linear')),
    c('sd_error', 'cor_error'))
  text = utils::capture.output(print(summary(names_fit(data))))
  expect_true('Structural parameters:' %in% text)
  expect_true(any(grepl('^  sd_re_selection +SD of the person effect in the selection equation$', text)))
  expect_true(any(grepl('^  cor_error +correlation of the selection error and the outcome error', text)))
})

test_that('the public argument lists are the documented ones', {
  expect_identical(names(formals(panel_select)), c('outcome', 'selection', 'data', 'id', 'family',
    'start', 'fixed', 'control', 'subset', 'na_action', 'init', 'vcov'))
  expect_identical(eval(formals(panel_select)$na_action), c('fail', 'omit', 'exclude'))
  expect_identical(eval(formals(panel_select)$init),
    c('default', 're', 'pooled', 'cross_sectional', 'zero', 'unif', 'norm'))
  expect_identical(eval(formals(panel_select)$vcov), c('opg', 'observed', 'sandwich'))
  expect_identical(names(formals(panel_control)), c('quadrature', 'nodes', 'quadrature_options',
    'observation_quadrature', 'observation_nodes', 'observation_options', 'mode_maxit', 'mode_tol',
    'exact_reductions', 'optimizer', 'method', 'maxit', 'reltol', 'factr', 'parscale', 'trace',
    'stationarity_warning', '...'))
  expect_identical(eval(formals(panel_control)$quadrature), c('gh', 'agh', 'gl'))
  expect_identical(eval(formals(panel_control)$observation_quadrature), c('gh', 'agh', 'gl'))
  expect_identical(eval(formals(predict.panelselect)$type),
    c('response', 'selection_probability', 'joint', 'selected_mean', 'link'))
  expect_identical(names(formals(check_quadrature)), c('object', 'nodes', 'observation_nodes',
    'tolerance', 'refit', 'compare', 'quadrature_options', '...'))
  expect_identical(eval(formals(check_quadrature)$compare), c('higher_order', 'gl'))
  expect_identical(names(formals(vcov.panelselect)), c('object', 'type', 'vcov', '...'))
})

test_that('staged checks and the observation time budget use the public names', {
  data = data.frame(id = rep(1:3, each = 2), d = c(1,1,1,0,0,0), y = c(1,1,0,NA,NA,NA))
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0,
    sd_re_selection = 1, sd_re_outcome = 1, cor_re = 0, cor_error = 0)
  staged = panel_select(y ~ 1, d ~ 1, data, 'id', 'probit', start = theta,
    fixed = theta[names(theta) != 'cor_error'],
    control = panel_control(quadrature = 'gl', nodes = 20, maxit = 0))
  expect_error(check_quadrature(staged, nodes = 20),
    'At least one effective person-effect order must increase and none may decrease', fixed = TRUE)
  counts = data.frame(id = rep(1:3, each = 2), count = c(0, 1, 2, 1, 0, 3))
  observation_only = panel_select(count ~ 1, NULL, counts, 'id', 'pln',
    fixed = c('outcome.(Intercept)' = 0, sd_re_outcome = 0, sd_error = 0.5),
    control = panel_control(observation_quadrature = 'gl', observation_nodes = 10))
  expect_error(check_quadrature(observation_only, nodes = 15),
    'This model has no person-effect integral to refine', fixed = TRUE)
  tick = 0
  budget = staged_inner_budget(list(group = c(1L, 1L), n_panels = 1L),
    list(seconds_per_panel = 5, seconds_total = 8), function() tick)
  budget$begin(1L); tick = 6
  expect_error(budget$check(), paste('PLN observation quadrature per-panel time budget exhausted',
    '(observation_options$seconds_per_panel = 5 s)'), fixed = TRUE)
  budget$end(check_limit = FALSE)
})

test_that('check rows name the refined integrals; saved checks are read in these names', {
  data = data.frame(id = rep(1:2, each = 3), y = c(0, 1, 2, 1, 0, 3))
  fit = panel_select(y ~ 1, NULL, data, 'id', 'pln',
    fixed = c('outcome.(Intercept)' = .2, sd_re_outcome = .3, sd_error = .2),
    control = panel_control(nodes = 3L, observation_nodes = 5L))
  check = check_quadrature(fit, nodes = 5L, observation_nodes = 9L)
  expect_identical(check$comparison$rule, c('person_effect', 'observation', 'both'))
  # Historical checks named the rows outer, inner and joint.
  old = check; old$comparison$rule = c('outer', 'inner', 'joint')
  expect_identical(current_check(old), check)
  expect_identical(utils::capture.output(print(old)), utils::capture.output(print(check)))
})
