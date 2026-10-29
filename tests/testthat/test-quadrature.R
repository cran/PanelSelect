test_that('normal rules integrate normalized moments including zero dimensions', {
  empty = normal_rule(3L, 0L)
  expect_equal(dim(empty$nodes), c(1L, 0L))
  expect_equal(empty$log_weights, 0)
  for (q in 1:2) {
    rule = normal_rule(c(5L, 7L), q)
    weights = exp(rule$log_weights)
    expect_equal(sum(weights), 1, tolerance = 1e-13)
    expect_equal(drop(crossprod(weights, rule$nodes^2)), rep(1, q), tolerance = 1e-13)
  }
})

test_that('adaptive scores differentiate recomputed modes and curvature', {
  for (model in core_models) for (nodes in c(1L, 3L)) {
    expect_core_score(core_fixture(model), core_control('agh', nodes = nodes, inner = 9L),
      tolerance = 5e-6)
  }
})

test_that('adaptive values approach independently evaluated high-order fixed rules', {
  expect_false(is.null(core_references))
  for (model in core_models) {
    fixture = core_fixture(model)
    adaptive = evaluate_model(fixture$theta, fixture$prepared,
      core_control('agh', nodes = 9L, inner = 25L))
    reference = core_reference(fixture, core_references$refined[[model]])
    expect_equal(fixture$theta, reference$theta)
    expect_lt(abs(adaptive$loglik-reference$loglik), 3e-6)
    expect_lt(max(abs(unname(adaptive$score)-unname(reference$score))), 3e-5)
  }
})

test_that('adaptive restricted models and empty-outcome panels remain differentiable', {
  for (model in core_models[startsWith(core_models, 'selection_')]) {
    fixture = core_fixture(model)
    names = grep('^cor_', names(fixture$theta), value = TRUE)
    fixture$theta[names] = 0
    expect_core_score(fixture, core_control('agh', nodes = 3L, inner = 9L), fixture$theta[names])
    fixture$prepared$d[] = 0; fixture$prepared$y[] = NA_real_
    expect_core_score(fixture, core_control('agh', nodes = 3L, inner = 9L))
  }
})

test_that('adaptive evaluation is repeatable across parameter call order', {
  fixture = core_fixture('selection_pln')
  control = core_control('agh', nodes = 3L, inner = 9L)
  first = evaluate_model(fixture$theta, fixture$prepared, control)
  changed = fixture$theta; changed['cor_re'] = 0.2
  evaluate_model(changed, fixture$prepared, control)
  second = evaluate_model(fixture$theta, fixture$prepared, control)
  expect_identical(second$loglik, first$loglik)
  expect_identical(second$score, first$score)
})

test_that('posterior weights and moving effect nodes have correct derivatives', {
  for (model in c('poisson_re', 'selection_gaussian', 'selection_pln')) {
    fixture = core_fixture(model)
    for (method in c('gh', 'agh')) {
      control = core_control(method, nodes = 3L, inner = 9L)
      posterior = posterior_rule(fixture$theta, fixture$prepared, control)[[1]]
      weights = exp(posterior$log_weights)
      expect_equal(sum(weights), 1, tolerance = 1e-12)
      # weight_score stores derivatives before posterior normalization.
      centered = sweep(posterior$weight_score, 2,
        drop(crossprod(weights, posterior$weight_score)), '-')
      expect_lt(max(abs(drop(crossprod(weights, centered)))), 1e-12)
      h = 1e-5
      for (j in seq_along(fixture$theta)) {
        plus = minus = fixture$theta; plus[j] = plus[j]+h; minus[j] = minus[j]-h
        upper = posterior_rule(plus, fixture$prepared, control)[[1]]
        lower = posterior_rule(minus, fixture$prepared, control)[[1]]
        expect_lt(max(abs(centered[, j]-(upper$log_weights-lower$log_weights)/(2*h))), 2e-5)
        expect_lt(max(abs(posterior$node_gradient[, , j]-
          (upper$nodes-lower$nodes)/(2*h))), 2e-5)
      }
    }
  }
})

test_that('nonfinite adaptive trial derivatives produce a recoverable integration error', {
  fixture = core_fixture('pln_re')
  control = core_control('agh', nodes = 3L, inner = 5L)
  before = evaluate_model(fixture$theta, fixture$prepared, control)
  for (parameter in c('outcome.(Intercept)', 'sd_re_outcome')) {
    trial = fixture$theta
    trial[parameter] = if (parameter == 'sd_re_outcome') 1e200 else 1000
    expect_error(evaluate_model(trial, fixture$prepared, control),
      'nonfinite likelihood derivatives', class = 'panelselect_integration_error')
  }
  after = evaluate_model(fixture$theta, fixture$prepared, control)
  expect_identical(after$loglik, before$loglik)
  expect_identical(after$score, before$score)
})
