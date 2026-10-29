test_that('conditional fixed selection preserves each tensor rule and score', {
  for (model in c('selection_probit', 'selection_poisson', 'selection_pln')) {
    fixture = core_fixture(model)
    fixed = if (model == 'selection_poisson') numeric() else c(cor_error = 0)
    fixture$theta[names(fixed)] = fixed
    for (rho in c(-.6, 0, .6)) for (nodes in list(5L, c(3L, 7L))) {
      fixture$theta['cor_re'] = rho
      control = core_control(nodes = nodes)
      generic = control; generic$exact_reductions = FALSE
      actual = evaluate_model(fixture$theta, fixture$prepared, control, fixed)
      reference = evaluate_model(fixture$theta, fixture$prepared, generic, fixed)
      expect_equal(actual, reference, tolerance = 1e-11)
      expect_equal(evaluate_model(fixture$theta, fixture$prepared, control, fixed,
        score = FALSE)$loglik, actual$loglik, tolerance = 1e-12)
      expect_core_score(fixture, control, fixed, tolerance = 3e-7)
      prepared = prepare_integration(fixture$prepared, control)
      settings = prepared$integration_cache$control
      context = make_model_context(fixture$theta, prepared, settings, fixed)
      for (batch in c(1L, 2L)) {
        blocked = integrate_selection_fixed(fixture$theta, prepared, context,
          settings$.rules$outer[[3L]], TRUE, batch = batch)
        blocked$score[, names(fixture$theta) %in% names(fixed)] = 0
        expect_equal(blocked$loglik, actual$loglik, tolerance = 1e-12)
        expect_equal(unname(blocked$score), unname(actual$score), tolerance = 1e-11)
      }
    }
  }
})

test_that('fixed probit reduction retains offsets, factors and selection-only panels', {
  data = data.frame(id = rep(1:6, each = 3), x = sin(1:18),
    f = factor(rep(c('a', 'b', 'c'), 6)), o = cos(1:18)/5,
    d = rep(c(0, 1, 1, 1, 0, 1), 3), y = rep(c(0, 1), 9))
  data$d[1:3] = 0; data$y[data$d == 0] = NA_real_
  prepared = prepare_data(y ~ x+f+offset(o), d ~ x+f+offset(-o),
    data, 'id', 'probit', 'fail')
  info = parameter_layout(prepared)
  theta = stats::setNames(rep(.1, nrow(info)), info$name)
  theta[c('sd_re_selection', 'sd_re_outcome', 'cor_re', 'cor_error')] = c(.5, .7, -.4, 0)
  fixed = c(cor_error = 0)
  control = core_control(nodes = c(4L, 7L))
  # Formula environments contain the test's changing local bindings. Compare
  # the numerical inputs, rather than serializing that enclosing environment.
  inputs = function() list(theta, control, fixed, prepared[c('x', 'w', 'y', 'd',
    'offset_outcome', 'offset_selection', 'group', 'group_ptr')])
  before = serialize(inputs(), NULL)
  generic = control; generic$exact_reductions = FALSE
  actual = evaluate_model(theta, prepared, control, fixed)
  expect_equal(actual, evaluate_model(theta, prepared, generic, fixed), tolerance = 1e-11)
  expect_equal(unname(actual$score[1L, !startsWith(names(theta), 'selection.') &
    names(theta) != 'sd_re_selection']), rep(0, sum(!startsWith(names(theta), 'selection.') &
      names(theta) != 'sd_re_selection')))
  expect_identical(serialize(inputs(), NULL), before)
})

test_that('generic paths remain available outside fixed two-effect tau-zero reduction', {
  fixture = core_fixture('selection_probit'); fixture$theta['cor_error'] = 0
  control = core_control(nodes = 3L)
  local_mocked_bindings(integrate_selection_fixed = function(...)
    stop('conditional fixed route reached'))
  expect_error(evaluate_model(fixture$theta, fixture$prepared, control,
    c(cor_error = 0)), 'conditional fixed route reached')
  # A free tau evaluated at zero retains its derivative and the joint kernel.
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control))
  expect_silent(posterior_rule(fixture$theta, fixture$prepared, control, c(cor_error = 0)))
  control$quadrature = 'agh'
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control, c(cor_error = 0)))
  control$quadrature = 'gh'
  fixture$theta['sd_re_outcome'] = 0; fixture$theta['cor_re'] = 0
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control,
    c(cor_error = 0, sd_re_outcome = 0, cor_re = 0)))
})

test_that('single-column grouped designs preserve ordinary arithmetic and empty groups', {
  x = cbind(c(1, Inf, 3), c(NA_real_, 2, 4))
  weight = matrix(c(0, 0, 2), 3L, 1L); ptr = c(0, 0, 1, 2, 3, 3)
  before = serialize(list(x, weight, ptr), NULL)
  expect_equal(group_weighted_design(x, weight, ptr),
    group_weighted_design_r(x, weight, ptr))
  expect_identical(serialize(list(x, weight, ptr), NULL), before)
})

test_that('zero GH edge weights retain their zero mass in the conditional rule', {
  fixture = core_fixture('selection_probit'); fixture$theta['cor_error'] = 0
  fixed = c(cor_error = 0); control = core_control(nodes = 3L)
  prepared = prepare_integration(fixture$prepared, control)
  context = make_model_context(fixture$theta, prepared, control, fixed)
  rule = normal_rule(3L, 2L)
  # Emulate underflow of the first marginal weight without a costly huge rule.
  weights = log(c(0, .8, .2))
  rule$log_weights = as.vector(outer(weights, weights, '+'))
  reduced = integrate_selection_fixed(fixture$theta, prepared, context, rule, TRUE)
  kernel = make_panel_kernel(fixture$theta, prepared, control, fixed,
    order = 1L, parameters = TRUE, latent = FALSE, context = context)
  full = kernel$evaluate(rule$nodes)
  reference = logsum_score_cpp(full$value+rule$log_weights, full$gradient, 1L)
  expect_equal(reduced$loglik, unname(reference[, 1L]), tolerance = 1e-12)
  expect_equal(reduced$score, reference[, -1L, drop = FALSE], tolerance = 1e-11)
})

# Retain the prior conditioned finite rule, including derivatives of its nodes.
linear_fixed_reference = function(theta, prepared, control, fixed, base, weights, score) {
  context = make_model_context(theta, prepared, control, fixed)
  context$marginal = NULL
  kernel = make_panel_kernel(theta, prepared, control, fixed, order = as.integer(score),
    parameters = score, latent = FALSE, context = context)
  n = prepared$n_panels
  result = kernel$evaluate(base[rep(seq_len(nrow(base)), each = n), , drop = FALSE])
  derivative = if (score) result$gradient[, kernel$parameter_columns, drop = FALSE] else
    matrix(0, n*nrow(base), length(theta))
  reduced = logsum_score_cpp(result$value+rep(weights, each = n), derivative, n)
  list(loglik = unname(reduced[, 1L]), score = unname(reduced[, -1L, drop = FALSE]))
}

test_that('Linear fixed moments preserve the conditioned rule and all first scores', {
  fixture = core_fixture('selection_gaussian'); fixture$theta['cor_error'] = 0
  control = core_control(nodes = 7L); fixed = c(cor_error = 0)
  rule = normal_rule(7L, 1L)
  for (rho in c(-.6, 0, .6)) {
    fixture$theta['cor_re'] = rho
    context = make_model_context(fixture$theta, fixture$prepared, control, fixed)
    for (score in c(FALSE, TRUE)) for (batch in c(1L, 3L, 7L)) {
      state = NULL
      for (first in seq.int(1L, 7L, by = batch)) {
        indices = seq.int(first, min(7L, first+batch-1L))
        base = rule$nodes[indices, , drop = FALSE]; weights = rule$log_weights[indices]
        actual = probit_re_linear_re_fixed_block(fixture$theta, fixture$prepared,
          context, base, weights, score)
        reference = linear_fixed_reference(fixture$theta, fixture$prepared,
          control, fixed, base, weights, score)
        expect_equal(actual, reference, tolerance = 1e-12, ignore_attr = TRUE)
        state = merge_quadrature(state, actual$loglik, actual$score)
      }
      full = probit_re_linear_re_fixed_block(fixture$theta, fixture$prepared,
        context, rule$nodes, rule$log_weights, score)
      expect_equal(state, full, tolerance = 1e-12)
    }
    expect_core_score(fixture, control, fixed, tolerance = 3e-7)
  }
})

test_that('Linear moments retain unequal panels, offsets, empty groups and zero weights', {
  data = data.frame(id = c(1, 1, 2, 3, 3, 3), x = sin(1:6), o = cos(1:6)/5,
    d = c(0, 0, 1, 1, 0, 1), y = c(NA, NA, .2, .4, NA, -.1))
  for (selection in list(d ~ x+offset(o), d ~ 0+offset(o))) {
    prepared = prepare_data(y ~ x+offset(-o), selection, data, 'id', 'gaussian', 'fail')
    info = parameter_layout(prepared)
    theta = stats::setNames(rep(.1, nrow(info)), info$name)
    theta[c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error')] = c(.5, .7, .4, .3, 0)
    control = core_control(nodes = 3L); fixed = c(cor_error = 0)
    context = make_model_context(theta, prepared, control, fixed)
    rule = normal_rule(3L, 1L); weights = log(c(0, .8, .2))
    before = serialize(list(theta, context, control, fixed), NULL)
    actual = probit_re_linear_re_fixed_block(theta, prepared, context, rule$nodes, weights, TRUE)
    expect_equal(actual, linear_fixed_reference(theta, prepared, control, fixed,
      rule$nodes, weights, TRUE), tolerance = 1e-12, ignore_attr = TRUE)
    state = NULL
    for (j in 1:3) {
      value = probit_re_linear_re_fixed_block(theta, prepared, context,
        rule$nodes[j, , drop = FALSE], weights[j], TRUE)
      state = merge_quadrature(state, value$loglik, value$score)
    }
    expect_equal(state, actual, tolerance = 1e-12)
    expect_identical(serialize(list(theta, context, control, fixed), NULL), before)
    fit_value = evaluate_model(theta, prepared, control, fixed)
    outcome = !startsWith(names(theta), 'selection.') & names(theta) != 'sd_re_selection'
    expect_true(all(fit_value$score[1L, outcome] == 0))
  }
})

test_that('Linear fixed moments do not replace adaptive, posterior or free-tau routes', {
  fixture = core_fixture('selection_gaussian'); fixture$theta['cor_error'] = 0
  control = core_control(nodes = 3L); fixed = c(cor_error = 0)
  local_mocked_bindings(probit_re_linear_re_fixed_block = function(...)
    stop('Linear fixed block reached'))
  expect_error(evaluate_model(fixture$theta, fixture$prepared, control, fixed), 'Linear fixed block reached')
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control))
  expect_silent(posterior_rule(fixture$theta, fixture$prepared, control, fixed))
  control$quadrature = 'agh'
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control, fixed))
  control$quadrature = 'gh'; control$exact_reductions = FALSE
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control, fixed))
  control$exact_reductions = TRUE; fixture$theta[c('sd_re_selection', 'cor_re')] = 0
  expect_silent(evaluate_model(fixture$theta, fixture$prepared, control,
    c(fixed, sd_re_selection = 0, cor_re = 0)))
})
