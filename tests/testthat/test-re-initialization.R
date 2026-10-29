# Family-specific auxiliary policies; numerical references are independent of
# the production Gaussian marginal and its optimizer-coordinate derivatives.
re_initial_fixture = function(family = 'pln', repeated = TRUE, selection = TRUE, cross = FALSE) {
  local_test_seed(2501)
  n = 48L
  data = data.frame(id = if (repeated) rep(1:12, each = 4) else seq_len(n),
    x = rnorm(n), z = rnorm(n), group = factor(rep(c('a', 'b', 'c'), 16)),
    off = seq(-.3, .2, length.out = n), selection_off = seq(.1, -.2, length.out = n),
    d = rep(c(1, 0, 1, 1, 0, 1, 1, 0), 6))
  effect = rep(rnorm(12, sd = .4), each = 4)
  data$y = if (family == 'gaussian') .2 + .3*data$x + effect + rnorm(n, sd = .6) else
    if (family == 'probit') as.numeric(.2 + .3*data$x + effect + rnorm(n) > 0) else
      rpois(n, exp(.2 + .2*data$x + effect))
  if (selection) data$y[data$d == 0] = NA_real_
  prepared = prepare_data(y ~ poly(x, 2) + group + offset(off),
    if (selection) d ~ z + offset(selection_off) else NULL,
    data, if (cross) NULL else 'id', family)
  list(data = data, prepared = prepared, info = parameter_layout(prepared))
}

re_initial_block = function(prepared, equation) {
  design = if (equation == 'selection') prepared$w else prepared$x
  c(paste0(equation, '.', colnames(design)), paste0('sd_re_', equation),
    if (equation == 'outcome' && prepared$family %in% c('gaussian', 'pln')) 'sd_error')
}

test_that('count RE auxiliaries preserve prepared designs, offsets and both integration levels', {
  for (family in c('poisson', 'pln')) {
    f = re_initial_fixture(family)
    supplied = c(selection.z = .25, 'outcome.poly(x, 2)1' = -.3)
    if (family == 'pln') supplied['sd_error'] = .4
    fixed = c('outcome.groupb' = -.17, sd_re_selection = .75)
    f$info$fixed = f$info$name %in% names(fixed)
    control = fit_control(panel_control(quadrature = 'agh', nodes = c(3, 5),
      observation_quadrature = 'gl', observation_nodes = 12,
      observation_options = list(reference_nodes = 17, bound = 9, max_nodes = 30000),
      optimizer = 'maxLik', method = 'NR', maxit = 0, reltol = 1e-4,
      mode_maxit = 17, mode_tol = 1e-7), init = 're', vcov = 'sandwich')
    control$integration_mesh = list(epoch = 99L)
    pooled = initial_parameters(f$prepared, f$info, supplied, fixed, fit_control(init = 'pooled'))
    calls = list()
    local_mocked_bindings(panel_select = function(outcome, selection, data, id, family,
        start, fixed, control, init, vcov, ...) {
      prepared = prepare_data(outcome, selection, data, id, family)
      info = parameter_layout(prepared)
      value = setNames(seq_len(nrow(info))/10 + .2, info$name)
      value[info$type == 'sd'] = .8
      value[names(fixed)] = fixed
      calls[[length(calls) + 1L]] <<- list(prepared = prepared, start = start, fixed = fixed,
        control = control, init = init, vcov = vcov, coefficients = value)
      # The mapping must use names, not the auxiliary vector's incidental order.
      list(coefficients = rev(value), optimizer = list(converged = TRUE))
    }, .package = 'PanelSelect')
    actual = initial_parameters(f$prepared, f$info, supplied, fixed, control)
    expect_length(calls, 2L)
    for (j in 1:2) {
      equation = c('selection', 'outcome')[j]
      rows = if (j == 1L) seq_len(f$prepared$n_rows) else which(f$prepared$selected)
      design = if (j == 1L) f$prepared$w else f$prepared$x
      response = if (j == 1L) f$prepared$d else f$prepared$y
      offset = if (j == 1L) f$prepared$offset_selection else f$prepared$offset_outcome
      block = re_initial_block(f$prepared, equation)
      aux = calls[[j]]
      expect_equal(unname(aux$prepared$x), unname(design[rows, , drop = FALSE]), tolerance = 0)
      expect_identical(aux$prepared$y, response[rows])
      expect_identical(aux$prepared$offset_outcome, offset[rows])
      expect_equal(aux$prepared$id_values[aux$prepared$group], f$prepared$group[rows])
      expect_false(aux$prepared$has_selection)
      expect_identical(aux$prepared$family, c('probit', family)[j])
      expect_identical(aux$control$quadrature, 'agh')
      expect_identical(aux$control$nodes, c(3L, 5L)[j])
      expect_identical(aux$control$mode_maxit, 17L)
      expect_identical(aux$control$mode_tol, 1e-7)
      expect_null(aux$control$integration_mesh)
      expect_identical(aux$control$optimizer, 'optim')
      expect_identical(aux$control$method, 'BFGS')
      expect_identical(aux$control$maxit, panel_control()$maxit)
      expect_identical(aux$control$reltol, panel_control()$reltol)
      expect_identical(aux$init, 'pooled')
      expect_identical(aux$vcov, 'opg')
      expect_equal(unname(aux$start), unname(pooled[block]))
      expect_equal(unname(aux$fixed), unname(fixed[intersect(block, names(fixed))]))
      expected = setNames(unname(aux$coefficients), block)
      expected[intersect(block, names(supplied))] = supplied[intersect(block, names(supplied))]
      expect_equal(actual[block], expected)
    }
    if (family == 'pln') {
      expect_identical(calls[[2L]]$control$observation_quadrature, 'gl')
      expect_identical(calls[[2L]]$control$observation_nodes, 12L)
      expect_identical(calls[[2L]]$control$observation_options, control$observation_options)
    }
    correlation = intersect(c('cor_re', 'cor_error'), names(actual))
    expect_identical(actual[correlation], pooled[correlation])
    expect_identical(actual[names(supplied)], supplied)
    expect_identical(actual[names(fixed)], fixed)
  }
})

test_that('invalid auxiliary results fall back only for the affected equation', {
  f = re_initial_fixture()
  supplied = c(selection.z = .29, 'outcome.poly(x, 2)1' = -.2, sd_error = .4)
  fixed = c('selection.(Intercept)' = .7, 'outcome.groupb' = .1)
  f$info$fixed = f$info$name %in% names(fixed)
  pooled = initial_parameters(f$prepared, f$info, supplied, fixed, fit_control(init = 'pooled'))
  for (failed in c('selection', 'outcome')) for (defect in c('error', 'nonconverged',
      'unnamed', 'missing', 'duplicate', 'unknown', 'nonfinite', 'negative_sd',
      'zero_free_sd', 'changed_fixed', 'character', 'matrix')) {
    calls = 0L
    local_mocked_bindings(panel_select = function(outcome, selection, data, id, family, fixed, ...) {
      calls <<- calls + 1L
      info = parameter_layout(prepare_data(outcome, NULL, data, id, family))
      value = setNames(rep(.46, nrow(info)), info$name)
      value[names(fixed)] = fixed
      converged = TRUE
      if (calls == match(failed, c('selection', 'outcome'))) {
        if (defect == 'error') stop('deliberate auxiliary error')
        if (defect == 'nonconverged') converged = FALSE
        if (defect == 'unnamed') value = unname(value)
        if (defect == 'missing') value = value[-1L]
        if (defect == 'duplicate') names(value)[2L] = names(value)[1L]
        if (defect == 'unknown') names(value)[1L] = 'unexpected_parameter'
        if (defect == 'nonfinite') value[1L] = Inf
        if (defect == 'negative_sd') value['sd_re_outcome'] = -.1
        if (defect == 'zero_free_sd') value['sd_re_outcome'] = 0
        if (defect == 'changed_fixed') value[names(fixed)[1L]] = fixed[[1L]] + .1
        if (defect == 'character') value = setNames(as.character(value), names(value))
        if (defect == 'matrix') value = matrix(value, nrow = 1L, dimnames = list(NULL, names(value)))
      }
      list(coefficients = value, optimizer = list(converged = converged))
    }, .package = 'PanelSelect')
    result = with_warnings(initial_parameters(f$prepared, f$info, supplied, fixed,
      fit_control(init = 're')))
    expect_identical(calls, 2L)
    expect_length(result$warnings, 1L)
    expect_match(result$warnings, paste(failed, 'RE initialization failed; using pooled starts:'))
    failed_block = re_initial_block(f$prepared, failed)
    expect_identical(result$value[failed_block], pooled[failed_block], info = defect)
    other = re_initial_block(f$prepared, if (failed == 'selection') 'outcome' else 'selection')
    other = setdiff(other, union(names(supplied), names(fixed)))
    expect_equal(unname(result$value[other]), rep(.46, length(other)), info = defect)
    expect_identical(result$value[names(supplied)], supplied)
    expect_identical(result$value[names(fixed)], fixed)
  }
})

test_that('converged count auxiliaries are accepted despite covariance and Q warnings', {
  f = re_initial_fixture()
  local_mocked_bindings(panel_select = function(outcome, selection, data, id, family, ...) {
    info = parameter_layout(prepare_data(outcome, NULL, data, id, family))
    warning('The requested covariance estimate is unavailable; fixture', call. = FALSE)
    warning('Optimizer stopped with Q_OPG = 100', call. = FALSE)
    list(coefficients = setNames(rep(.6, nrow(info)), info$name),
      optimizer = list(converged = TRUE, stationarity = list(Q_OPG = 100, warning = TRUE)),
      covariance_status = list(opg = list(valid = FALSE)))
  }, .package = 'PanelSelect')
  expect_silent({ actual = initial_parameters(f$prepared, f$info, NULL, numeric(), fit_control(init = 're')) })
  parameters = union(re_initial_block(f$prepared, 'selection'), re_initial_block(f$prepared, 'outcome'))
  expect_equal(unname(actual[parameters]), rep(.6, length(parameters)))
  expect_identical(unname(actual[c('cor_re', 'cor_error')]), c(0, 0))
})

test_that('the Gaussian auxiliary receives original designs, offsets and natural restrictions', {
  f = re_initial_fixture('gaussian')
  pooled = initial_parameters(f$prepared, f$info, NULL, numeric(), fit_control(init = 'pooled'))
  supplied = c(pooled[re_initial_block(f$prepared, 'selection')],
    'outcome.poly(x, 2)1' = .23, sd_re_outcome = .65)
  fixed = c('outcome.groupb' = -.2, sd_error = .7)
  f$info$fixed = f$info$name %in% names(fixed)
  pooled = initial_parameters(f$prepared, f$info, supplied, fixed, fit_control(init = 'pooled'))
  control = fit_control(panel_control(quadrature = 'agh', nodes = c(4, 7),
    optimizer = 'maxLik', method = 'NR', maxit = 0, reltol = 1e-4), init = 're')
  call = NULL
  local_mocked_bindings(panel_select = function(...) stop('complete selection block must be skipped'),
    linear_re_initial_fit = function(prepared, start, fixed, control) {
      call <<- list(prepared = prepared, start = start, fixed = fixed, control = control)
      value = setNames(rep(.8, length(start)), names(start)); value[names(fixed)] = fixed
      list(coefficients = rev(value), optimizer = list(converged = TRUE))
    }, .package = 'PanelSelect')
  actual = initial_parameters(f$prepared, f$info, supplied, fixed, control)
  block = re_initial_block(f$prepared, 'outcome')
  expect_identical(call$prepared, f$prepared)
  expect_identical(call$start, pooled[block])
  expect_identical(call$fixed, fixed)
  expect_identical(call$control$nodes, 7L)
  expect_identical(call$control$quadrature, 'agh')
  expect_identical(call$control$optimizer, 'optim')
  expect_identical(call$control$method, 'BFGS')
  expect_identical(call$control$maxit, panel_control()$maxit)
  expect_identical(call$control$reltol, panel_control()$reltol)
  expect_identical(actual[names(supplied)], supplied)
  expect_identical(actual[names(fixed)], fixed)
  free = setdiff(block, union(names(supplied), names(fixed)))
  expect_equal(unname(actual[free]), rep(.8, length(free)))
})

test_that('singleton and zero-scale skips follow the auxiliary family identification', {
  cases = list(
    list(family = 'probit', fixed = numeric(), expected = character()),
    list(family = 'probit', fixed = c(sd_re_outcome = .6), expected = 'probit'),
    list(family = 'poisson', fixed = numeric(), expected = 'poisson'),
    list(family = 'poisson', fixed = c(sd_re_outcome = 0), expected = character()),
    list(family = 'gaussian', fixed = numeric(), expected = character()),
    list(family = 'gaussian', fixed = c(sd_re_outcome = 0), expected = 'gaussian'),
    list(family = 'gaussian', fixed = c(sd_error = .5), expected = 'gaussian'),
    list(family = 'pln', fixed = numeric(), expected = character()),
    list(family = 'pln', fixed = c(sd_re_outcome = 0), expected = 'pln'),
    list(family = 'pln', fixed = c(sd_error = 0), expected = 'pln'),
    list(family = 'pln', fixed = c(sd_re_outcome = 0, sd_error = 0), expected = character()))
  calls = character()
  local_mocked_bindings(panel_select = function(outcome, selection, data, id, family, fixed, ...) {
    calls <<- c(calls, family)
    info = parameter_layout(prepare_data(outcome, NULL, data, id, family))
    value = setNames(rep(.6, nrow(info)), info$name); value[names(fixed)] = fixed
    list(coefficients = value, optimizer = list(converged = TRUE))
  }, linear_re_initial_fit = function(prepared, start, fixed, control) {
    calls <<- c(calls, 'gaussian')
    list(coefficients = start, optimizer = list(converged = TRUE))
  }, .package = 'PanelSelect')
  for (case in cases) {
    f = re_initial_fixture(case$family, repeated = FALSE)
    fixed = resolve_fixed(c(sd_re_selection = 0, case$fixed), f$info, f$prepared)
    f$info$fixed = f$info$name %in% names(fixed)
    calls = character()
    expect_silent({ actual = initial_parameters(f$prepared, f$info, NULL, fixed, fit_control(init = 're')) })
    expect_identical(calls, case$expected, info = case$family)
    expect_identical(actual[names(fixed)], fixed)
  }
  # A fixed positive selection SD permits its otherwise unidentified singleton probit fit.
  f = re_initial_fixture('probit', repeated = FALSE)
  fixed = resolve_fixed(c(sd_re_selection = .7, sd_re_outcome = 0), f$info, f$prepared)
  f$info$fixed = f$info$name %in% names(fixed)
  calls = character()
  initial_parameters(f$prepared, f$info, NULL, fixed, fit_control(init = 're'))
  expect_identical(calls, 'probit')
})

test_that('standalone, cross-sectional and complete-start routes never recurse into auxiliaries', {
  local_mocked_bindings(panel_select = function(...) stop('unexpected standalone auxiliary'),
    linear_re_initial_fit = function(...) stop('unexpected Gaussian auxiliary'), .package = 'PanelSelect')
  fixtures = c(lapply(c('probit', 'poisson', 'pln'), re_initial_fixture, selection = FALSE),
    lapply(c('probit', 'gaussian', 'pln'), re_initial_fixture, cross = TRUE))
  for (f in fixtures) {
    pooled = initial_parameters(f$prepared, f$info, NULL, numeric(), fit_control(init = 'pooled'))
    for (method in c('default', 're', 'pooled', 'cross_sectional'))
      expect_identical(initial_parameters(f$prepared, f$info, NULL, numeric(),
        fit_control(init = method)), pooled)
  }
  for (family in c('probit', 'gaussian', 'poisson', 'pln')) {
    f = re_initial_fixture(family)
    supplied = initial_parameters(f$prepared, f$info, NULL, numeric(), fit_control(init = 'pooled'))
    expect_identical(initial_parameters(f$prepared, f$info, supplied, numeric(), fit_control(init = 're')), supplied)
  }
})

# Selected panels of sizes three, two, one and zero, with an outcome offset.
re_gaussian_reference_data = function() {
  data = data.frame(id = rep(1:8, each = 4), x = rep(c(-.8, .1, .7, 1.2), 8),
    off = seq(-.2, .2, length.out = 32), d = rep(c(1, 1, 0, 1), 8))
  data$d[c(5L, 9L, 10L, 29:32)] = 0
  person = rep(c(-1.2, -.8, -.4, -.1, .2, .5, .9, 1.1), each = 4)
  error = rep(c(.2, -.3, .1, 0), 8)
  data$y = 1 + .5*data$x + data$off + person + error
  data$y[data$d == 0] = NA_real_
  prepare_data(y ~ x + offset(off), d ~ 1, data, 'id', 'gaussian')
}

re_gaussian_loglik_reference = function(theta, prepared) {
  vapply(seq_len(prepared$n_panels), function(i) {
    rows = which(prepared$group == i & prepared$selected)
    n = length(rows)
    if (!n) return(0)
    covariance = diag(theta[['sd_error']]^2, n) +
      matrix(theta[['sd_re_outcome']]^2, n, n)
    mean = drop(prepared$x[rows, , drop = FALSE] %*% theta[paste0('outcome.', colnames(prepared$x))]) +
      prepared$offset_outcome[rows]
    residual = prepared$y[rows] - mean
    -.5*(n*log(2*pi) + as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
      sum(residual*solve(covariance, residual)))
  }, numeric(1))
}

test_that('Gaussian auxiliary likelihood and scores match a dense independent covariance', {
  prepared = re_gaussian_reference_data()
  theta = c('outcome.(Intercept)' = .4, outcome.x = -.2, sd_re_outcome = .7, sd_error = .5)
  actual = linear_outcome_marginal(theta, prepared, conditional = FALSE)
  expect_equal(actual$value, re_gaussian_loglik_reference(theta, prepared), tolerance = 1e-12)
  for (name in names(theta)) {
    plus = minus = theta; h = 1e-6
    plus[name] = plus[name] + h; minus[name] = minus[name] - h
    numerical = (re_gaussian_loglik_reference(plus, prepared) -
      re_gaussian_loglik_reference(minus, prepared))/(2*h)
    expect_equal(actual$gradient[, match(name, names(theta))], numerical, tolerance = 1e-7, info = name)
  }
  expect_identical(actual$value[8L], 0)
  expect_identical(unname(actual$gradient[8L, ]), rep(0, length(theta)))
})

test_that('Gaussian auxiliary optimization agrees with the independent exact likelihood', {
  prepared = re_gaussian_reference_data()
  start = c('outcome.(Intercept)' = .4, outcome.x = .2, sd_re_outcome = .8, sd_error = .6)
  natural = function(eta) { value = eta; value[3:4] = exp(eta[3:4]); value }
  eta = start; eta[3:4] = log(eta[3:4])
  reference = stats::optim(eta, function(value) -sum(re_gaussian_loglik_reference(natural(value), prepared)),
    method = 'BFGS', control = list(reltol = 1e-12, maxit = 1000))
  local_mocked_bindings(information_matrix = function(...) stop('auxiliary requested information'),
    integrate_prepared = function(...) stop('exact Gaussian auxiliary requested quadrature'), .package = 'PanelSelect')
  fit = linear_re_initial_fit(prepared, start,
    control = panel_control(reltol = 1e-12, maxit = 1000, parscale = 1))
  expect_identical(reference$convergence, 0L)
  expect_true(fit$optimizer$converged)
  expect_equal(fit$loglik, -reference$value, tolerance = 1e-9)
  expect_equal(fit$coefficients, natural(reference$par), tolerance = 1e-5)
  expect_null(fit$optimizer$information)
  expect_null(fit$covariance)
  fixed = c(outcome.x = .4, sd_re_outcome = 0, sd_error = .6)
  start[names(fixed)] = fixed
  fit = linear_re_initial_fit(prepared, start, fixed,
    panel_control(reltol = 1e-12, parscale = 1))
  rows = which(prepared$selected)
  intercept = mean(prepared$y[rows] - prepared$offset_outcome[rows] - .4*prepared$x[rows, 'x'])
  expect_true(fit$optimizer$converged)
  expect_equal(unname(fit$coefficients['outcome.(Intercept)']), intercept, tolerance = 1e-8)
  expect_identical(fit$coefficients[names(fixed)], fixed)
})

test_that('initialization methods require an exact scalar name in fit and update calls', {
  data = data.frame(id = 1:6, y = c(0, 1, 2, 1, 3, 2))
  fit = panel_select(y ~ 1, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0),
    init = 'pooled', control = panel_control(maxit = 0))
  invalid = list('ProbitRE', 'def', 'cross', 'cross-sectional', 'r', 'pool', 'z', 'RE', ' re', 're ', NULL, NA_character_,
    character(), c('re', 'pooled'), eval(formals(panel_select)$init), factor('re'), 1,
    matrix('re', 1L, 1L))
  for (method in invalid) {
    expect_error(panel_select(y ~ 1, NULL, data, 'id', 'poisson',
      fixed = c(sd_re_outcome = 0), init = method, control = panel_control(maxit = 0)), 'init')
    expect_error(update(fit, init = method), 'init')
    expect_error(update(fit, init = method, evaluate = FALSE), 'init')
  }
  for (method in c('default', 'pooled', 'zero', 'unif', 'norm'))
    expect_identical(update(fit, init = method, evaluate = FALSE)$init, method)
})

test_that('saved ProbitRE fits display but require an explicit supported initializer to update', {
  data = data.frame(id = 1:6, y = c(0, 1, 2, 1, 3, 2))
  formula = stats::as.formula('y ~ 1', env = baseenv())
  fit = panel_select(formula, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0),
    init = 'pooled', control = panel_control(maxit = 0))
  previous = fit; previous$control$init = 'ProbitRE'
  stored = serialize(previous, NULL)
  expect_output(print(previous), 'PanelSelect')
  expect_silent(invisible(summary(previous)))
  expect_identical(coef(previous), coef(fit))
  expect_identical(serialize(previous, NULL), stored)
  expect_error(update(previous), 'init')
  expect_error(update(previous, evaluate = FALSE), 'init')
  expect_identical(update(previous, init = 're', evaluate = FALSE)$init, 'default')
  updated = update(previous, init = 're')
  expect_identical(updated$control$init, 'default')
  expect_identical(coef(updated), coef(fit))
})
