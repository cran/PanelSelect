split_single_fixture = function() {
  data = data.frame(id = rep(1L, 4), d = c(0, 1, 0, 0),
    y = c(NA, 0, NA, NA), x = c(-.83, .27, .91, 1.31),
    os = c(.043, -.031, .071, -.057), oy = c(0, .051, 0, 0))
  prepared = prepare_data(y ~ x + offset(oy), d ~ x + offset(os), data, 'id', 'probit', 'fail')
  theta = c('selection.(Intercept)' = .173, selection.x = -.283,
    'outcome.(Intercept)' = -.237, outcome.x = .193,
    sd_re_selection = 1.37, sd_re_outcome = 3.13, cor_re = -.31, cor_error = 0)
  list(prepared = prepared, theta = theta, panel = prepare_split_panels(prepared)$panels[[1]])
}

# Independently construct the NEW finite conditional rule, without production
# cut/node/weight/score helpers. The midpoint map is algebraically equivalent.
split_single_finite_reference = function(theta, panel, order = 3L, bound = 5) {
  a = drop(panel$w %*% theta[paste0('selection.', colnames(panel$w))]) + panel$offset_selection
  b = drop(panel$x %*% theta[paste0('outcome.', colnames(panel$x))]) + panel$offset_outcome
  delta = theta[['sd_re_selection']]; lambda = theta[['sd_re_outcome']]; rho = theta[['cor_re']]
  denominator = sqrt(1 + lambda^2*(1 - rho^2))
  h = b/denominator; g = lambda*rho/denominator
  levels = c(-4, -2, -1, 0, 1, 2, 4)
  cuts = c(-bound, seq(-bound, bound, by = 2), bound,
    as.vector(outer(-a, levels, '+')/delta))
  cuts = sort(unique(cuts[is.finite(cuts) & abs(cuts) <= bound]))
  rule = statmod::gauss.quad(order, kind = 'legendre')
  probability = vapply(seq_len(length(cuts) - 1L), function(j) {
    width = cuts[j + 1L] - cuts[j]
    u = (cuts[j] + cuts[j + 1L])/2 + width/2*rule$nodes
    selection = vapply(u, function(z) prod(pnorm(panel$selection_sign*(a + delta*z))), numeric(1))
    sum(width/2*rule$weights*dnorm(u)*selection*pnorm(panel$outcome_sign*(h + g*u)))
  }, numeric(1))
  log(sum(probability))
}

test_that('single-outcome analytic scores differentiate its complete moving finite rule', {
  z = split_single_fixture(); options = list(bound = 5, max_nodes = 1000000L)
  for (rho in c(-.31, 0, .27)) {
    theta = z$theta; theta['cor_re'] = rho
    actual = split_single_outcome_panel(theta, z$panel, 3L, options, fixed = c(cor_error = 0))
    expect_equal(actual$loglik, split_single_finite_reference(theta, z$panel), tolerance = 1e-13)
    expect_identical(actual$loglik, split_single_outcome_panel(theta, z$panel, 3L, options,
      FALSE, fixed = c(cor_error = 0))$loglik)
    for (step in c(1e-5, 2.5e-6)) {
      finite = theta*0
      for (name in setdiff(names(theta), 'cor_error')) {
        h = step*(1 + abs(theta[name])); plus = minus = theta
        plus[name] = plus[name] + h; minus[name] = minus[name] - h
        finite[name] = (split_single_finite_reference(plus, z$panel) -
          split_single_finite_reference(minus, z$panel))/(2*h)
      }
      expect_lt(max(abs(actual$score - finite)/(1 + abs(finite))), 2e-6)
    }
  }
})

test_that('single-outcome values and zero-loading rho score match independent integrals', {
  z = split_single_fixture(); theta = z$theta
  reference = function(theta, derivative = FALSE) {
    a = drop(z$panel$w %*% theta[1:2]) + z$panel$offset_selection
    b = drop(z$panel$x %*% theta[3:4]) + z$panel$offset_outcome
    denominator = sqrt(1 + theta['sd_re_outcome']^2*(1 - theta['cor_re']^2))
    integrand = function(u) vapply(u, function(node) {
      selection = prod(pnorm(z$panel$selection_sign*(a + theta['sd_re_selection']*node)))
      argument = z$panel$outcome_sign*(b + theta['sd_re_outcome']*theta['cor_re']*node)/denominator
      if (derivative) dnorm(node)*selection*dnorm(argument)*z$panel$outcome_sign*
        theta['sd_re_outcome']*node/denominator else dnorm(node)*selection*pnorm(argument)
    }, numeric(1))
    value = integrate(integrand, -Inf, Inf, rel.tol = 1e-11, abs.tol = 1e-13, subdivisions = 1000)
    expect_identical(value$message, 'OK')
    value$value
  }
  options = list(bound = 12, max_nodes = 1000000L)
  for (lambda in c(.8, 8, 30)) for (rho in c(-.95, -.3, -1e-8, 0, 1e-8, .3, .95)) {
    theta[c('sd_re_outcome', 'cor_re')] = c(lambda, rho)
    actual = split_single_outcome_panel(theta, z$panel, 20L, options, fixed = c(cor_error = 0))
    expect_equal(actual$loglik, log(reference(theta)), tolerance = 1e-9)
    if (rho == 0) {
      exact_score = reference(theta, TRUE)/reference(theta)
      expect_equal(unname(actual$score['cor_re']), unname(exact_score), tolerance = 1e-9)
      expect_gt(abs(actual$score['cor_re']), 1e-4)
    }
    expect_identical(actual$dimension, 1L)
  }
})

test_that('single-outcome zero-scale restrictions retain exact constrained reductions', {
  z = split_single_fixture(); options = list(bound = 12, max_nodes = 1000000L)
  control = panel_control(quadrature = 'gl', nodes = 20)
  for (name in c('sd_re_selection', 'sd_re_outcome')) {
    theta = z$theta; theta[name] = 0
    fixed = c(cor_error = 0, setNames(0, name))
    actual = evaluate_split_binary(theta, z$prepared, control, fixed)
    theta['cor_re'] = -theta['cor_re']
    other = evaluate_split_binary(theta, z$prepared, control, fixed)
    expect_identical(actual$loglik, other$loglik)
    expect_equal(actual$score[1, c(name, 'cor_re')], c(0, 0), ignore_attr = TRUE)
    expect_error(split_single_outcome_panel(theta, z$panel, 20L, options,
      fixed = c(cor_error = 0)), 'explicitly fixed')
    if (name == 'sd_re_selection') {
      a = drop(z$panel$w %*% theta[1:2]) + z$panel$offset_selection
      b = drop(z$panel$x %*% theta[3:4]) + z$panel$offset_outcome
      exact = sum(pnorm(z$panel$selection_sign*a, log.p = TRUE)) +
        pnorm(z$panel$outcome_sign*b/sqrt(1 + theta['sd_re_outcome']^2), log.p = TRUE)
      expect_equal(actual$loglik, unname(exact), tolerance = 1e-13)
    }
  }
})

test_that('transition-level zero-loading limits keep one continuous finite objective', {
  z = split_single_fixture(); panel = z$panel; panel$offset_outcome[] = 0
  options = list(bound = 5, max_nodes = 1000000L)
  theta = z$theta; theta['cor_re'] = theta['outcome.x'] = 0
  lambda = theta[['sd_re_outcome']]; d0 = sqrt(1 + lambda^2)
  for (level in c(-2, -1, 0, 1, 2)) {
    theta['outcome.(Intercept)'] = level*d0
    base = split_single_outcome_panel(theta, panel, 3L, options, fixed = c(cor_error = 0))
    for (k in c(-.731, .419)) {
      path = function(t) {
        at = theta; at['cor_re'] = t
        d = sqrt(1 + lambda^2*(1 - t^2)); g = lambda*t/d
        at['outcome.(Intercept)'] = d*(level + k*g)
        at
      }
      # b'(0)=k*lambda; selection geometry has no singular outcome cuts.
      derivative = base$score['cor_re'] + k*lambda*base$score['outcome.(Intercept)']
      for (step in c(1e-5, 2.5e-6)) {
        right = split_single_outcome_panel(path(step), panel, 3L, options, FALSE,
          fixed = c(cor_error = 0))$loglik
        left = split_single_outcome_panel(path(-step), panel, 3L, options, FALSE,
          fixed = c(cor_error = 0))$loglik
        expect_equal(unname(derivative), (right - left)/(2*step), tolerance = 2e-6)
      }
      for (t in c(-1e-8, 1e-8)) expect_lt(abs(split_single_outcome_panel(path(t), panel,
        3L, options, FALSE, fixed = c(cor_error = 0))$loglik - base$loglik), 1e-6)
    }
  }
})

test_that('single-outcome route is explicit, preserves fixed masks and keeps all native gates', {
  z = split_single_fixture(); options = list(bound = 12, max_nodes = 1000000L)
  expect_error(split_single_outcome_panel(z$theta, z$panel, 12L, options), 'explicitly fixed')
  expect_error(split_physical_panel(z$theta, z$panel, 12L, options), 'explicitly fixed')
  expect_error(split_single_outcome_panel(z$theta, z$panel, 12L,
    list(bound = 12, max_nodes = 1L), fixed = c(cor_error = 0)), 'node budget')
  control = panel_control(quadrature = 'gl', nodes = 12, quadrature_options = list(reference_nodes = 15))
  actual = evaluate_split_binary(z$theta, z$prepared, control, c(cor_error = 0))
  expect_identical(actual$diagnostics$routes, 'conditional_single_outcome')
  expect_equal(actual$score[, 'cor_error'], 0, ignore_attr = TRUE)
  original = evaluate_split_binary(z$theta, z$prepared, control, c(cor_error = 0), single_outcome = FALSE)
  expect_identical(original$diagnostics$routes, 'physical')
  expect_lt(actual$diagnostics$nodes, original$diagnostics$nodes)
  expect_equal(actual$loglik, original$loglik, tolerance = 1e-9)
  expect_error(evaluate_split_binary(z$theta, z$prepared,
    panel_control(quadrature = 'gl', nodes = 20, quadrature_options = list(bound = 3)),
    c(cor_error = 0), FALSE), 'tail bound', class = 'panelselect_split_accuracy_error')
  coarse = panel_control(quadrature = 'gl', nodes = 3, quadrature_options = list(reference_nodes = 2,
    panel_loglik_tol = 1, total_loglik_tol = 1, panel_score_tol = 1e-12, total_score_tol = 1e-12))
  expect_error(evaluate_split_binary(z$theta, z$prepared, coarse, c(cor_error = 0)),
    'score refinement|branch ambiguity', class = 'panelselect_split_accuracy_error')
})

test_that('conditional cut crossings carry a branch bound with unequal coefficient designs', {
  z = split_single_fixture(); theta = z$theta; theta[1:4] = 0
  panel = z$panel; panel$offset_selection[] = panel$offset_outcome[] = 0
  options = list(bound = 8, max_nodes = 1000000L)
  coarse = split_single_outcome_panel(theta, panel, 2L, options, fixed = c(cor_error = 0))
  expect_gt(max(coarse$ambiguity), 1e-4)
  for (name in c('selection.x', 'outcome.x')) {
    step = 1e-7; plus = minus = theta; plus[name] = step; minus[name] = -step
    right = (split_single_outcome_panel(plus, panel, 2L, options, FALSE,
      fixed = c(cor_error = 0))$loglik - coarse$loglik)/step
    left = (coarse$loglik - split_single_outcome_panel(minus, panel, 2L, options, FALSE,
      fixed = c(cor_error = 0))$loglik)/step
    expect_lte(max(abs(c(left, right) - coarse$score[name])), unname(coarse$ambiguity[name]) + 1e-5)
  }
})
