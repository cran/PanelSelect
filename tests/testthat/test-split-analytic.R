split_analytic_fixture = function() {
  data = data.frame(id = rep(1L, 3), d = c(1, 0, 1), y = c(0, NA, 1),
    x = c(-.67, .23, 1.17), os = c(.043, -.031, .071), oy = c(.051, 0, -.021))
  prepared = prepare_data(y ~ x + offset(oy), d ~ x + offset(os), data, 'id', 'probit', 'fail')
  theta = c('selection.(Intercept)' = .173, selection.x = -.283,
    'outcome.(Intercept)' = -.237, outcome.x = .193,
    sd_re_selection = 1.37, sd_re_outcome = 2.13, cor_re = -.31, cor_error = 0)
  list(prepared = prepared, theta = theta,
    panel = prepare_split_panels(prepared)$panels[[1]])
}

# Independent low-order finite sum: no production node, weight, derivative or
# density helper is used. Finite differences are regression tests only.
split_analytic_reference = function(theta, panel, order = 2L, bound = 5) {
  axis = function(index, scale) {
    cuts = sort(unique(c(-bound, seq(-bound, bound, by = 2), bound,
      as.vector(outer(-index, c(-4, -2, -1, 0, 1, 2, 4), '+')/scale))))
    cuts = cuts[abs(cuts) <= bound]
    rule = statmod::gauss.quad(order, kind = 'legendre')
    list(nodes = unlist(lapply(seq_len(length(cuts) - 1L), function(j)
      (cuts[j] + cuts[j + 1L])/2 + diff(cuts[c(j, j + 1L)])/2*rule$nodes)),
      weights = unlist(lapply(diff(cuts)/2, function(width) width*rule$weights)))
  }
  a = drop(panel$w %*% theta[paste0('selection.', colnames(panel$w))]) + panel$offset_selection
  b = drop(panel$x %*% theta[paste0('outcome.', colnames(panel$x))]) + panel$offset_outcome
  u = axis(a, theta['sd_re_selection']); v = axis(b, theta['sd_re_outcome'])
  product = function(nodes, index, scale, signs)
    apply(outer(nodes, index, function(z, a) z*scale + a), 1L,
      function(z) prod(pnorm(signs*z)))
  s = product(u$nodes, a, theta['sd_re_selection'], panel$selection_sign)
  y = product(v$nodes, b, theta['sd_re_outcome'], panel$outcome_sign)
  rho = theta['cor_re']; q = 1 - rho^2
  density = outer(u$nodes, v$nodes, function(u, v)
    exp(-(u^2 + v^2 - 2*rho*u*v)/(2*q))/(2*pi*sqrt(q)))
  log(sum(density*outer(u$weights*s, v$weights*y)))
}

test_that('split analytical scores differentiate the actual moving finite rule', {
  z = split_analytic_fixture()
  options = list(bound = 5, max_nodes = 1000000L)
  actual = split_physical_panel(z$theta, z$panel, 2L, options, fixed = c(cor_error = 0))
  expect_equal(actual$loglik, split_analytic_reference(z$theta, z$panel), tolerance = 1e-13)
  expect_identical(actual$loglik,
    split_physical_panel(z$theta, z$panel, 2L, options, score = FALSE, fixed = c(cor_error = 0))$loglik)
  for (step in c(1e-5, 2.5e-6)) {
    finite = z$theta*0
    for (name in setdiff(names(z$theta), 'cor_error')) {
      h = step*(1 + abs(z$theta[name])); plus = minus = z$theta
      plus[name] = plus[name] + h; minus[name] = minus[name] - h
      finite[name] = (split_analytic_reference(plus, z$panel) -
        split_analytic_reference(minus, z$panel))/(2*h)
    }
    expect_lt(max(abs(actual$score - finite)/(1 + abs(finite))), 1e-6)
  }
})

test_that('a free zero rho retains its cross-moment score', {
  z = split_analytic_fixture(); z$theta['cor_re'] = 0
  factor = function(index, signs, scale) {
    f = function(node) vapply(node, function(u) dnorm(u)*prod(pnorm(signs*(index + scale*u))), numeric(1))
    probability = integrate(f, -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-12)$value
    mean = integrate(function(u) u*f(u), -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-12)$value/probability
    c(probability = probability, mean = mean)
  }
  a = drop(z$panel$w %*% z$theta[1:2]) + z$panel$offset_selection
  b = drop(z$panel$x %*% z$theta[3:4]) + z$panel$offset_outcome
  s = factor(a, z$panel$selection_sign, z$theta['sd_re_selection'])
  y = factor(b, z$panel$outcome_sign, z$theta['sd_re_outcome'])
  actual = evaluate_split_binary(z$theta, z$prepared,
    panel_control(quadrature = 'gl', nodes = 20), c(cor_error = 0))
  expect_equal(actual$loglik, unname(log(s['probability']*y['probability'])), tolerance = 1e-9)
  expect_equal(unname(actual$score[1, 'cor_re']), unname(s['mean']*y['mean']), tolerance = 1e-9)
  expect_gt(abs(actual$score[1, 'cor_re']), 1e-5)
  expect_equal(actual$score[1, 'cor_error'], 0, ignore_attr = TRUE)
})

test_that('exact cut crossings have bounded branch ambiguity and an explicit gate', {
  z = split_analytic_fixture(); z$theta[1:4] = 0
  z$panel$offset_selection[] = z$panel$offset_outcome[] = 0
  z$prepared$offset_selection[] = z$prepared$offset_outcome[] = 0
  options = list(bound = 8, max_nodes = 1000000L)
  coarse = split_physical_panel(z$theta, z$panel, 2L, options, fixed = c(cor_error = 0))
  expect_gt(max(coarse$ambiguity), 1e-3)
  for (name in c('selection.x', 'outcome.x')) {
    h = 1e-7; plus = minus = z$theta; plus[name] = h; minus[name] = -h
    right = (split_physical_panel(plus, z$panel, 2L, options, FALSE,
      fixed = c(cor_error = 0))$loglik - coarse$loglik)/h
    left = (coarse$loglik - split_physical_panel(minus, z$panel, 2L, options, FALSE,
      fixed = c(cor_error = 0))$loglik)/h
    expect_lte(max(abs(c(left, right) - coarse$score[name])), unname(coarse$ambiguity[name]) + 1e-5)
  }
  control = panel_control(quadrature = 'gl', nodes = 3,
    quadrature_options = list(reference_nodes = 2, panel_loglik_tol = 1, total_loglik_tol = 1,
      panel_score_tol = 1e-10, total_score_tol = 1e-10))
  expect_error(evaluate_split_binary(z$theta, z$prepared, control, c(cor_error = 0)),
    'analytical score|branch ambiguity', class = 'panelselect_split_accuracy_error')
  fine = evaluate_split_binary(z$theta, z$prepared,
    panel_control(quadrature = 'gl', nodes = 20), c(cor_error = 0))
  expect_true(fine$diagnostics$score_checked)
  expect_identical(fine$diagnostics$gradient_method, 'analytic')
})

test_that('fixed zero scales and selection-only panels reduce without phantom outcomes', {
  z = split_analytic_fixture(); control = panel_control(quadrature = 'gl', nodes = 20)
  for (scale in c('sd_re_selection', 'sd_re_outcome')) {
    theta = z$theta; theta[scale] = 0; theta['cor_re'] = .8
    fixed = c(cor_error = 0, setNames(0, scale))
    actual = evaluate_split_binary(theta, z$prepared, control, fixed)
    theta['cor_re'] = -.8
    expect_identical(actual$loglik, evaluate_split_binary(theta, z$prepared, control, fixed)$loglik)
    expect_equal(actual$score[1, c(scale, 'cor_re')], c(0, 0), ignore_attr = TRUE)
    expect_error(evaluate_split_binary(theta, z$prepared, control, c(cor_error = 0)), 'explicitly fixed')
  }
  none = z$prepared; none$d[] = 0; none$y[] = NA_real_
  actual = evaluate_split_binary(z$theta, none, control, c(cor_error = 0))
  other = z$theta; other[c('outcome.(Intercept)', 'outcome.x', 'sd_re_outcome', 'cor_re')] = c(20, -30, 17, -.98)
  expect_identical(actual$loglik, evaluate_split_binary(other, none, control, c(cor_error = 0))$loglik)
  expect_equal(actual$score[1, c(names(other)[3:4], 'sd_re_outcome', 'cor_re')], rep(0, 4), ignore_attr = TRUE)
})

test_that('split caching ignores structural missing outcomes and preserves canonical values', {
  z = split_analytic_fixture(); control = panel_control(quadrature = 'gl', nodes = 20)
  altered = z$prepared; missing = which(altered$d == 0)
  altered$y[missing] = Inf; altered$x[missing, ] = NA_real_; altered$offset_outcome[missing] = NA_real_
  expect_identical(prepare_split_panels(altered), prepare_split_panels(z$prepared))
  actual = evaluate_split_binary(z$theta, z$prepared, control, c(cor_error = 0))
  altered$split_cache = prepare_split_panels(altered)
  cached = evaluate_split_binary(z$theta, altered, control, c(cor_error = 0))
  expect_identical(actual$loglik, cached$loglik)
  expect_identical(actual$score, cached$score)
  value = evaluate_split_binary(z$theta, altered, control, c(cor_error = 0), score = FALSE)
  expect_identical(value$loglik, actual$loglik)
  expect_null(value$score)
  expect_false(value$diagnostics$score_checked)
  expect_identical(value$diagnostics$check_scope, 'likelihood/refinement/tail only')
  expect_error(evaluate_split_binary(z$theta[-1], z$prepared, control, c(cor_error = 0)), 'do not match')
})
