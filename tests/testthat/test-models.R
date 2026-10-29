test_that('all seven fixed rules agree with independent likelihoods and scores', {
  expect_false(is.null(core_references))
  for (model in core_models) for (sign in c(-1, 1)) {
    fixture = core_fixture(model, sign)
    control = core_control(exact_reductions = FALSE)
    actual = evaluate_model(fixture$theta, fixture$prepared, control)
    reference = core_reference(fixture, core_references$fixed[[paste(model, sign)]])
    expect_equal(fixture$theta, reference$theta)
    expect_equal(fixture$panel, reference$panel)
    expect_equal(actual$loglik, reference$loglik, tolerance = 2e-10, info = model)
    expect_lt(max(abs(unname(actual$score)-unname(reference$score))), 2e-9)
  }
})

test_that('fixed scores differentiate the actual finite rule', {
  for (model in core_models) expect_core_score(core_fixture(model),
    core_control(nodes = 3L, inner = 5L), tolerance = 3e-6)
})

test_that('restricted likelihoods retain correct free scores', {
  for (model in core_models[startsWith(core_models, 'selection_')]) {
    for (restriction in list('cor_re', 'cor_error', c('cor_re', 'cor_error'))) {
      fixture = core_fixture(model)
      restriction = intersect(restriction, names(fixture$theta))
      if (!length(restriction)) next
      fixture$theta[restriction] = 0
      fixed = fixture$theta[restriction]
      expect_core_score(fixture, core_control(nodes = 3L, inner = 5L), fixed)
      # Only Gaussian conditioning changes the finite grid. Other restrictions
      # factor the same grid exactly and therefore need no high-order rule.
      nodes = if (model == 'selection_gaussian') 19L else 5L
      control = core_control(nodes = nodes, inner = 9L)
      reduced = evaluate_model(fixture$theta, fixture$prepared, control, fixed)
      control$exact_reductions = FALSE
      general = evaluate_model(fixture$theta, fixture$prepared, control, fixed)
      expect_lt(max(abs(reduced$loglik-general$loglik)), 2e-6)
      free = !names(fixture$theta) %in% restriction
      expect_lt(max(abs(reduced$score[, free]-general$score[, free])), 3e-5)
    }
  }
})

test_that('a freely estimated zero correlation keeps its derivative', {
  for (model in core_models[startsWith(core_models, 'selection_')]) {
    fixture = core_fixture(model)
    correlations = grepl('^cor_', names(fixture$theta))
    fixture$theta[correlations] = 0
    result = expect_core_score(fixture, core_control(nodes = 5L, inner = 7L))
    expect_gt(max(abs(result$score[, correlations])), 1e-4)
  }
})

test_that('all-unselected panels have no outcome information', {
  for (model in core_models[startsWith(core_models, 'selection_')]) {
    fixture = core_fixture(model)
    fixture$prepared$d[] = 0
    fixture$prepared$selected[] = FALSE
    fixture$prepared$y[] = NA_real_
    for (method in c('gh', 'agh')) {
      control = core_control(method, nodes = 5L)
      result = evaluate_model(fixture$theta, fixture$prepared, control)
      outcome = grepl('^outcome|sd_re_outcome|sd_error|cor_', names(fixture$theta))
      expect_equal(unname(result$score[, outcome]), rep(0, sum(outcome)))
      fixture$prepared$y[] = 999
      other = evaluate_model(fixture$theta, fixture$prepared, control)
      expect_equal(other$loglik, result$loglik)
      expect_equal(other$score, result$score)
    }
  }
})

test_that('zero scale models use their ordinary likelihood without quadrature', {
  for (model in c('probit_re', 'poisson_re', 'pln_re')) {
    fixture = core_fixture(model)
    scales = grep('^sd_', names(fixture$theta), value = TRUE)
    fixture$theta[scales] = 0
    fixed = fixture$theta[scales]
    result = evaluate_model(fixture$theta, fixture$prepared, core_control('agh'), fixed)
    eta = drop(fixture$prepared$x %*% fixture$theta[1:2])
    expected = if (model == 'probit_re') sum(pnorm((2*fixture$prepared$y-1)*eta, log.p = TRUE)) else
      sum(dpois(fixture$prepared$y, exp(eta), log = TRUE))
    expect_equal(result$loglik, expected, tolerance = 1e-12)
    expect_equal(result$diagnostics$dimension, 0L)
  }
})

test_that('unbalanced batches preserve independent panel values and scores', {
  for (model in core_models) {
    fixture = core_fixture(model)
    first = fixture$prepared
    second = first
    for (name in c('x', 'w')) second[[name]] = first[[name]][c(3, 1), , drop = FALSE]
    for (name in c('y', 'd', 'selected', 'offset_outcome', 'offset_selection'))
      second[[name]] = first[[name]][c(3, 1)]
    if (second$has_selection) {
      second$d[] = 0; second$selected[] = FALSE; second$y[] = NA_real_
    }
    second$group = c(1L, 1L); second$group_ptr = c(0L, 2L); second$n_rows = 2L
    batch = first
    for (name in c('x', 'w')) batch[[name]] = rbind(first[[name]], second[[name]])
    for (name in c('y', 'd', 'selected', 'offset_outcome', 'offset_selection'))
      batch[[name]] = c(first[[name]], second[[name]])
    batch$group = c(1L, 1L, 1L, 2L, 2L); batch$group_ptr = c(0L, 3L, 5L)
    batch$n_rows = 5L; batch$n_panels = 2L; batch$id_values = c('first', 'second')
    for (method in c('gh', 'agh')) {
      control = core_control(method, nodes = 3L, inner = 5L)
      a = evaluate_model(fixture$theta, first, control)
      b = evaluate_model(fixture$theta, second, control)
      together = evaluate_model(fixture$theta, batch, control)
      expect_equal(together$loglik, c(a$loglik, b$loglik), tolerance = 1e-10)
      expect_equal(together$score, rbind(a$score, b$score), tolerance = 1e-10)
    }
  }
})
