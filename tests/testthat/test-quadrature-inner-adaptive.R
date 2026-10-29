inner_option_fixture = function(model = 'pln_re', large = FALSE) {
  selected = model != 'pln_re'; panel = model != 'probit_pln'
  theta = c('outcome.(Intercept)' = if (large) -.5 else .25)
  if (selected) theta = c('selection.(Intercept)' = if (large) 1.7 else -.3, theta)
  if (panel) theta = c(theta, sd_re_outcome = if (large) 2.5 else .7)
  if (selected && panel) theta = c(theta, sd_re_selection = if (large) 2.2 else .6,
    cor_re = if (large) .07 else .35)
  theta = c(theta, sd_error = if (large) 2.5 else .5)
  if (selected) theta = c(theta, cor_error = if (large) .4 else -.25)
  one = matrix(1, 1, 1, dimnames = list(NULL, '(Intercept)'))
  prepared = list(model_id = model, structure = if (panel) 'panel' else 'cross_sectional',
    family = 'pln', has_selection = selected, x = one,
    w = if (selected) one else matrix(numeric(), 1, 0), y = 10, d = 1,
    offset_outcome = 0, offset_selection = 0, group = 1L,
    group_ptr = c(0L, 1L), id_values = 'one', n_rows = 1L, n_panels = 1L, selected = TRUE)
  list(theta = theta, prepared = prepared)
}

# Independent singleton collapse: combine person and occasion log-rate effects,
# then condition the normal selection index on that single combined normal.
inner_option_reference = function(fixture) {
  theta = fixture$theta
  value = function(name) if (name %in% names(theta)) unname(theta[name]) else 0
  a = value('selection.(Intercept)'); b = value('outcome.(Intercept)')
  delta = value('sd_re_selection'); lambda = value('sd_re_outcome')
  gamma = value('sd_error'); rho = value('cor_re'); tau = value('cor_error')
  y = fixture$prepared$y
  stopifnot(length(y) == 1L)
  scale = sqrt(lambda^2 + gamma^2); covariance = delta*lambda*rho + gamma*tau
  kernel = function(z) {
    eta = b + scale*z
    out = y*eta-exp(eta)-lgamma(y+1)+dnorm(z, log = TRUE)
    if (fixture$prepared$has_selection) out = out + pnorm(
      (a + covariance*z/scale)/sqrt(1+delta^2-covariance^2/scale^2), log.p = TRUE)
    out
  }
  mode = optimize(function(z) -kernel(z), c(-14, 14))$minimum
  peak = kernel(mode)
  peak + log(integrate(function(z) exp(kernel(z)-peak), -14, 14,
    rel.tol = 1e-11, subdivisions = 500L)$value)
}

test_that('inner quadrature is explicit and the frozen default references are preserved', {
  expect_identical(panel_control()$observation_quadrature, 'gh')
  expect_error(panel_control(observation_quadrature = 'unknown'))
  # Setting only the observation order keeps the fixed GH observation rule.
  expect_identical(panel_control(quadrature = 'gh', nodes = 7L, observation_nodes = 13L)$observation_quadrature, 'gh')
  for (model in c('pln_re', 'selection_pln')) for (method in c('gh', 'agh')) {
    fixture = core_fixture(model)
    implicit = core_control(method)
    explicit = implicit; explicit$observation_quadrature = 'gh'
    first = evaluate_model(fixture$theta, fixture$prepared, implicit)
    second = evaluate_model(fixture$theta, fixture$prepared, explicit)
    expect_identical(first$loglik, second$loglik)
    expect_identical(first$score, second$score)
    # These references predate this option and are never regenerated here.
    # This fixture has fixed H7/I13 references, but no adaptive H7/I13 set.
    if (method == 'gh') {
      saved = core_reference(fixture, core_references$fixed[[paste(model, 1)]])
      expect_equal(sum(first$loglik), saved$loglik, tolerance = 2e-11)
      expect_equal(unname(colSums(first$score)), unname(saved$score), tolerance = 2e-9)
    }
  }
})

test_that('all three PLN routes expose independently checked adaptive inner values', {
  for (model in c('pln_re', 'probit_re_pln_re', 'probit_pln')) for (large in c(FALSE, TRUE)) {
    fixture = inner_option_fixture(model, large)
    control = panel_control(quadrature = 'agh', nodes = 20L,
      observation_nodes = 20L, observation_quadrature = 'agh')
    actual = evaluate_model(fixture$theta, fixture$prepared, control)
    reference = inner_option_reference(fixture)
    # Selected large-scale outer H20 retains ~4.6e-6 error; do not hide it.
    expect_lt(abs(sum(actual$loglik)-reference), if (large) 1e-5 else 1e-8)
    expect_identical(actual$loglik, evaluate_model(fixture$theta, fixture$prepared,
      control, score = FALSE)$loglik)
  }
})

test_that('adaptive inner quadrature resolves the motivating rare-count singleton', {
  fixture = inner_option_fixture('pln_re')
  fixture$prepared$y = 100
  control = panel_control(quadrature = 'agh', nodes = 20L, observation_nodes = 20L,
    observation_quadrature = 'agh')
  actual = evaluate_model(fixture$theta, fixture$prepared, control)
  reference = inner_option_reference(fixture)
  # Independent collapse combines both normal effects before numerical
  # integration; it does not reuse either package quadrature or its score.
  expect_lt(abs(sum(actual$loglik)-reference), 1e-8)
  expect_true(all(is.finite(actual$score)))
  expect_identical(actual$loglik, evaluate_model(fixture$theta, fixture$prepared,
    control, score = FALSE)$loglik)
})

test_that('complete nested scores and observed information differentiate the new objective', {
  for (model in c('pln_re', 'probit_re_pln_re', 'probit_pln')) {
    fixture = inner_option_fixture(model)
    info = parameter_layout(fixture$prepared)
    theta = fixture$theta[info$name]
    free = c('outcome.(Intercept)', 'sd_error', if (model != 'pln_re') 'cor_error')
    info$fixed = !info$name %in% free
    fixed = theta[info$fixed]
    for (method in c('gh', 'agh')) {
      control = panel_control(quadrature = method, nodes = 7L, observation_nodes = 13L,
        observation_quadrature = 'agh', mode_tol = 1e-11)
      eta = parameter_to_optimizer(theta, info)
      evaluator = fit_evaluator(theta, info, fixture$prepared, control, fixed)
      value = function(x) -sum(evaluator$evaluate(x, trial = FALSE)$loglik)
      actual = evaluator$evaluate(eta, trial = FALSE)
      numerical = core_difference(value, eta)
      expect_lt(max(abs(numerical + colSums(actual$optimizer_score))), 3e-6)
      information = information_matrix(eta, evaluator)
      # Independent second differences of the objective for every diagonal entry
      # and the cross term of the last two free coordinates.
      h = 1e-4; center = value(eta)
      cells = rbind(cbind(seq_along(eta), seq_along(eta)), c(length(eta)-1L, length(eta)))
      independent = numeric(nrow(cells))
      for (k in seq_len(nrow(cells))) {
        i = cells[k, 1L]; j = cells[k, 2L]
        ei = ej = rep(0, length(eta)); ei[i] = h; ej[j] = h
        independent[k] = if (i == j) (value(eta+ei)-2*center+value(eta-ei))/h^2 else
          (value(eta+ei+ej)-value(eta+ei-ej)-value(eta-ei+ej)+value(eta-ei-ej))/(4*h^2)
      }
      expect_lt(max(abs(information[cells]-independent)/(1+abs(independent))), 3e-4)
    }
  }
})

test_that('adaptive inner values compose with the outer rule on a two-row panel', {
  # Independent composition: the test-only adaptive PLN engine
  # (reference/pln-inner.R) integrates each row's observation error at every
  # node of a seven-node outer GH rule; the panel likelihood is their weighted sum.
  theta = c('outcome.(Intercept)' = .25, outcome.x = -.4, sd_re_outcome = .7, sd_error = .5)
  x = c(-.3, .8); y = c(3, 0)
  prepared = list(model_id = 'pln_re', structure = 'panel', family = 'pln',
    has_selection = FALSE, x = cbind('(Intercept)' = 1, x = x), w = matrix(numeric(), 2, 0),
    y = y, d = c(1, 1), offset_outcome = c(0, 0), offset_selection = c(0, 0),
    group = c(1L, 1L), group_ptr = c(0L, 2L), id_values = 'one', n_rows = 2L,
    n_panels = 1L, selected = c(TRUE, TRUE))
  control = panel_control(quadrature = 'gh', nodes = 7L, observation_nodes = 13L,
    observation_quadrature = 'agh', mode_tol = 1e-12)
  actual = evaluate_model(theta, prepared, control, score = FALSE)$loglik
  outer = ref_normal_rule(7L); inner = ref_normal_rule(13L)
  rule = list(nodes = matrix(inner$nodes, ncol = 1L), log_weights = log(inner$weights))
  row_logs = vapply(outer$nodes, function(z) {
    b = theta[['outcome.(Intercept)']] + theta[['outcome.x']]*x + theta[['sd_re_outcome']]*z
    sum(local_pln_adaptive(b, y, theta[['sd_error']], 0L, lgamma(y + 1), rule,
      mode_tol = 1e-12)$value)
  }, numeric(1))
  expect_equal(actual, log(sum(outer$weights*exp(row_logs))), tolerance = 1e-11)
})

test_that('public supplied starts and fixed controls retain the inner option through update', {
  data = data.frame(id = c(1,1,2,2), d = c(1,0,1,1), y = c(1,NA,2,0))
  for (model in c('pln_re', 'probit_re_pln_re', 'probit_pln')) {
    fixture = inner_option_fixture(model)
    input = data
    if (model == 'pln_re') input$y[is.na(input$y)] = 0
    control = panel_control(nodes = 3L, observation_nodes = 9L, observation_quadrature = 'agh')
    fit = panel_select(y ~ 1, if (model == 'pln_re') NULL else d ~ 1, input, if (model == 'probit_pln') NULL else 'id',
      'pln',
      start = fixture$theta, fixed = fixture$theta, control = control)
    expect_equal(coef(fit)[names(fixture$theta)], fixture$theta)
    expect_identical(fit$control$observation_quadrature, 'agh')
    updated = update(fit)
    expect_identical(updated$control$observation_quadrature, 'agh')
    expect_equal(as.numeric(logLik(updated)), as.numeric(logLik(fit)), tolerance = 1e-12)
    checked = check_quadrature(fit, nodes = 4L, observation_nodes = 11L)
    expect_true(all(checked$comparison$observation_quadrature == 'agh'))
    expect_null(checked$refit)
  }
})

test_that('inner boundary restrictions and numerical failures remain explicit', {
  rule = normal_rule(7L, 1L); rule$quadrature = 'agh'
  empty = integrate_pln_inner(NULL, numeric(), numeric(), .5, 0, 3L, numeric(),
    rule, pln_re_inner_kernel)
  expect_length(empty$value, 0L)
  expect_error(integrate_pln_inner(NULL, 1e6, 0, .5, 0, 1L, 0,
    rule, pln_re_inner_kernel), class = 'panelselect_integration_error')
  fixture = inner_option_fixture()
  info = parameter_layout(fixture$prepared); theta = fixture$theta[info$name]
  control = panel_control(observation_quadrature = 'agh', nodes = 3L)
  evaluator = fit_evaluator(theta, info, fixture$prepared, control, numeric())
  trial = parameter_to_optimizer(theta, info); trial['outcome.(Intercept)'] = 1e6
  expect_false(evaluator$evaluate(trial)$valid)
  # Repeated panels pass the earlier singleton-identification guard, so these
  # calls exercise their intended starting-value restrictions.
  data = data.frame(id = rep(1:3, each = 2), y = c(0,1,2,0,1,2))
  expect_error(panel_select(y ~ 1, NULL, data, 'id', 'pln',
    start = c(sd_error = 0), control = control), 'free standard deviations')
  expect_error(panel_select(y ~ 1, NULL, data, 'id', 'pln', start = c(sd_error = .5),
    fixed = c(sd_error = .6), control = control), 'conflicts with a fixed')
})

test_that('fast and higher inner derivatives agree across negative Mills tail cutoffs', {
  for (target in c(-25.001,-25,-24.999,-15.001,-15,-14.999,-4.001,-4,-3.999)) {
    # tau=0 makes the selection index equal a at every inner node; tau remains
    # a derivative coordinate, so its mixed derivatives are still checked.
    parameter = c(target,.25,.5,0)
    fast = pln_inner_adaptive_first(parameter, 2, 13L, TRUE)
    high = pln_inner_adaptive_derivatives(parameter, 2, 13L, TRUE)
    expect_equal(fast$value, high$value, tolerance = 2e-12)
    expect_lt(max(abs(fast$gradient-high$gradient)/(1+abs(high$gradient))), 2e-9)
  }
})
