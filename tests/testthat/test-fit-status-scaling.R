test_that('maxLik BFGS iteration limits remain nonconverged for an exact probit', {
  skip_if_not_installed('maxLik')
  data = data.frame(id = rep(1:10, each = 3), y = rep(c(1, 1, 0), 10))
  fit = NULL
  expect_warning({ fit = panel_select(y ~ 1, NULL, data, 'id', 'probit',
    start = c('outcome.(Intercept)' = 0), fixed = c(sd_re_outcome = 0),
    control = panel_control(quadrature = 'gh', optimizer = 'maxLik',
      method = 'BFGS', maxit = 1L, stationarity_warning = NULL), vcov = 'opg') },
    'Optimizer did not report convergence')
  expect_equal(fit$optimizer$code, 1)
  expect_false(fit$optimizer$converged)
  expect_false(fit$optimizer$information_computed)
  beta = unname(coef(fit)['outcome.(Intercept)'])
  exact = sum(pnorm((2*data$y - 1)*beta, log.p = TRUE))
  expect_equal(fit$loglik, exact, tolerance = 1e-12)
  expect_true(is.finite(fit$optimizer$gradient))
  expect_warning({ not_run = panel_select(y ~ 1, NULL, data, 'id', 'probit',
    start = c('outcome.(Intercept)' = 0), fixed = c(sd_re_outcome = 0),
    control = panel_control(quadrature = 'gh', optimizer = 'maxLik',
      method = 'BFGS', maxit = 0L, stationarity_warning = NULL), vcov = 'opg') }, NA)
  expect_identical(not_run$optimizer$status, 'not_run')
  expect_equal(not_run$optimizer$code, 0)
  expect_false(not_run$optimizer$converged)
  expect_equal(not_run$loglik, nrow(data)*log(.5), tolerance = 1e-12)
})

test_that('backend exit codes retain their method-specific meanings', {
  for (backend in c('optim', 'maxLik')) {
    methods = if (backend == 'optim') c('BFGS', 'CG', 'Nelder-Mead', 'L-BFGS-B') else
      c('BFGS', 'CG', 'NM', 'Nelder-Mead')
    for (method in methods) {
      control = panel_control(optimizer = backend, method = method)
      options = optimizer_options(control)
      expect_identical(optimizer_status(0, control, options, 1L), 'converged')
      for (code in c(1, 2, 8))
        expect_identical(optimizer_status(code, control, options, 1L), 'not_converged')
    }
    # ?optim: SANN always returns code 0, so code 0 is not evidence of convergence.
    control = panel_control(optimizer = backend, method = 'SANN')
    options = optimizer_options(control)
    for (code in c(0, 1, 2, 8))
      expect_identical(optimizer_status(code, control, options, 1L), 'no_convergence_test')
  }
  for (method in c('NR', 'Newton-Raphson', 'BHHH', 'BFGSR', 'BFGS-R')) {
    control = panel_control(optimizer = 'maxLik', method = method)
    options = optimizer_options(control)
    for (code in c(1, 2, 8))
      expect_identical(optimizer_status(code, control, options, 1L), 'converged')
    for (code in c(0, 3, 4))
      expect_identical(optimizer_status(code, control, options, 1L), 'not_converged')
  }
})

status_scaling_fixture = function() {
  data = data.frame(id = rep(1:12, each = 3),
    x = rep(c(-.7, .2, .9), 12) + rep(seq(-.4, .5, length.out = 12), each = 3),
    d = rep(c(1, 0, 1, 1, 1, 0), 6), y = rep(c(0, NA, 1, 1, 0, NA), 6))
  prepared = prepare_data(y ~ x, d ~ x, data, 'id', 'probit', 'fail')
  theta = c('selection.(Intercept)' = .17, selection.x = -.28,
    'outcome.(Intercept)' = -.24, outcome.x = .19,
    sd_re_selection = 1.37, sd_re_outcome = 2.13, cor_re = -.31, cor_error = 0)
  fixed = c(cor_error = 0)
  info = parameter_layout(prepared); info$fixed = info$name == 'cor_error'
  control = fit_control(panel_control(quadrature = 'gl', nodes = 20), vcov = 'opg')
  value = evaluate_model(theta, prepared, control, fixed)
  free = !info$fixed
  score = sweep(value$score[, free, drop = FALSE], 2L, parameter_jacobian(theta, info)[free], '*')
  list(data = data, prepared = prepared, theta = theta, fixed = fixed,
    info = info, control = control, value = value, score = score)
}

test_that('BFGS uses one initial free-panel-score metric and respects supplied scales', {
  z = status_scaling_fixture()
  # Coefficient scales follow their panel-score norms; log-SD and
  # atanh-correlation scales keep the cap at one.
  norm = sqrt(colSums(z$score^2))
  coefficient = grepl('^(selection|outcome)[.]', colnames(z$score))
  expected = pmax(sqrt(.Machine$double.eps), ifelse(coefficient, 1/norm, 1/pmax(1, norm)))
  for (backend in c('optim', 'maxLik')) {
    control = z$control; control$optimizer = backend
    automatic = optimizer_options(control, z$prepared, z$info, z$score)
    expect_equal(unname(automatic$parscale), unname(expected), tolerance = 1e-14)
    expect_identical(automatic$reltol, control$reltol)
    explicit = panel_control(quadrature = 'gl', nodes = 20,
      optimizer = backend, parscale = 2, reltol = 1e-10)
    actual = optimizer_options(explicit, z$prepared, z$info, z$score)
    expect_equal(unname(rep_len(actual$parscale, ncol(z$score))), rep(2, ncol(z$score)))
    expect_identical(actual$reltol, 1e-10)
    requested = setNames(seq_len(ncol(z$score))/3, rev(colnames(z$score)))
    explicit$parscale = requested
    actual = optimizer_options(validate_panel_control(explicit), z$prepared, z$info, z$score)
    expect_equal(unname(actual$parscale), unname(requested[colnames(z$score)]))
    expect_identical(explicit$parscale, requested)
  }
  explicit = z$control; explicit$parscale = 1
  expect_equal(unname(rep_len(optimizer_options(explicit, z$prepared, z$info,
    z$score)$parscale, ncol(z$score))), rep(1, ncol(z$score)))
  other = z$control; other$method = 'Nelder-Mead'
  expect_null(optimizer_options(other, z$prepared, z$info, z$score)$parscale)
  for (invalid in list(0, -1, NA_real_, Inf, 'automatic'))
    expect_error(panel_control(parscale = invalid), 'parscale')
  bad = z$control; bad$parscale = c(wrong = 2, other = 3)
  expect_error(optimizer_options(bad, z$prepared, z$info, z$score), 'parscale')
})

for (backend in c('optim', 'maxLik')) {
  test_that(paste('BFGS applies parscale to an exact objective without latent integration:', backend), {
    if (backend == 'maxLik') skip_if_not_installed('maxLik')
    y = 0:4; eta = c(beta = log(.8)); scale = .1
    evaluator = list(evaluate = function(value) {
      mu = exp(unname(value))
      list(loglik = dpois(y, mu, log = TRUE),
        optimizer_score = matrix(y - mu, length(y), 1L, dimnames = list(NULL, 'beta')))
    })
    # BFGS starts from an identity metric in scaled coordinates, so its first
    # accepted natural-coordinate step is parscale^2 times the total score.
    expected = eta + scale^2*sum(y - .8)
    control = fit_control(panel_control(optimizer = backend, method = 'BFGS',
      maxit = 1L, parscale = scale, stationarity_warning = NULL))
    result = optimize_model(eta, evaluator, control)
    expect_equal(result$par, expected, tolerance = 1e-12)
    expect_equal(result$convergence, 1)
    if (backend == 'maxLik') for (method in c('NR', 'BHHH', 'BFGSR'))
      expect_error(optimize_model(eta, evaluator, panel_control(optimizer = 'maxLik',
        method = method, parscale = scale)), 'parscale.*native maxLik')
  })
}

test_that('fixed scaling preserves the complete finite-rule objective and gradient chain', {
  z = status_scaling_fixture()
  eta = parameter_to_optimizer(z$theta, z$info)
  scaled = z$control; scaled$parscale = setNames(seq_along(eta)/4, names(eta))
  evaluator = fit_evaluator(z$theta, z$info, z$prepared, scaled, z$fixed)
  actual = evaluator$evaluate(eta, trial = FALSE)
  # Compare scaling at the identical optimizer-coordinate point. exp(log(sd))
  # and tanh(atanh(rho)) can differ from the original natural doubles by one ulp.
  unscaled = fit_evaluator(z$theta, z$info, z$prepared, z$control, z$fixed)$evaluate(eta, trial = FALSE)
  expect_identical(actual$loglik, unscaled$loglik)
  expect_identical(actual$score, unscaled$score)
  expect_identical(actual$optimizer_score, unscaled$optimizer_score)
  expect_identical(actual$diagnostics$canonical_order, 20L)
  expect_identical(actual$diagnostics$reference_order, 10L)
  # Differentiate the same finite likelihood independently of its score path.
  for (name in names(eta)) {
    plus = minus = eta; h = 2e-6
    plus[name] = plus[name] + h; minus[name] = minus[name] - h
    likelihood = function(e) sum(evaluate_model(parameter_from_optimizer(e,
      z$theta, z$info), z$prepared, scaled, z$fixed, score = FALSE)$loglik)
    numerical = (likelihood(plus) - likelihood(minus))/(2*h)
    expect_equal(unname(sum(actual$optimizer_score[, name])), numerical, tolerance = 2e-6)
  }
  # A fitted control is the validated numerical control with its fit settings.
  expect_identical(serialize(scaled, NULL),
    serialize(fit_control(scaled, scaled$na_action, scaled$init, scaled$vcov), NULL))
})

for (backend in c('optim', 'maxLik')) {
  test_that(paste('fixed scales leave reported coordinates and exact OPG covariance unchanged:', backend), {
    if (backend == 'maxLik') skip_if_not_installed('maxLik')
    data = data.frame(id = rep(1:6, each = 2),
      x = c(-1, 0, 1, 2, -.5, .5, -1, 1, 0, 2, -.5, .5),
      y = c(0, 1, 3, 1, 0, 2, 2, 0, 1, 4, 0, 1))
    theta = c('outcome.(Intercept)' = log(1.4), outcome.x = .15)
    design = cbind('(Intercept)' = 1, x = data$x)
    mu = exp(drop(design %*% theta))
    score = rowsum(design*as.numeric(data$y - mu), data$id, reorder = FALSE)
    covariance = solve(crossprod(score))
    for (scale in list(NULL, 1, c(outcome.x = .3, 'outcome.(Intercept)' = 2))) {
      fit = panel_select(y ~ x, NULL, data, 'id', 'poisson', start = theta,
        fixed = c(sd_re_outcome = 0), control = panel_control(optimizer = backend,
          maxit = 0, parscale = scale, stationarity_warning = NULL), vcov = 'opg')
      expect_equal(fit$loglik, sum(dpois(data$y, mu, log = TRUE)), tolerance = 1e-12)
      expect_identical(fit$optimizer$par, theta)
      expect_equal(unname(fit$score[, names(theta), drop = FALSE]), unname(score), tolerance = 1e-12)
      expect_equal(unname(fit$optimizer$gradient), unname(colSums(score)), tolerance = 1e-12)
      expect_equal(unname(vcov(fit)[names(theta), names(theta)]), unname(covariance), tolerance = 1e-12)
      expect_false(fit$optimizer$information_computed)
      expect_identical(fit$optimizer$status, 'not_run')
      expect_identical(fit$control$parscale, scale)
      expect_null(fit$optimizer$rejected)
    }
  })
}

for (backend in c('optim', 'maxLik')) {
  test_that(paste('BFGS reaches the exact Poisson maximum with fixed scales:', backend), {
    if (backend == 'maxLik') skip_if_not_installed('maxLik')
    data = data.frame(id = 1:5, y = 0:4)
    maximum = sum(dpois(data$y, 2, log = TRUE))
    for (scale in c(.1, 1)) {
      fit = panel_select(y ~ 1, NULL, data, 'id', 'poisson',
        start = c('outcome.(Intercept)' = log(.8)), fixed = c(sd_re_outcome = 0),
        control = panel_control(optimizer = backend, method = 'BFGS',
          parscale = scale, reltol = 1e-12, stationarity_warning = NULL), vcov = 'opg')
      expect_true(fit$optimizer$converged)
      expect_equal(unname(coef(fit)['outcome.(Intercept)']), log(2), tolerance = 1e-5)
      expect_equal(fit$loglik, maximum, tolerance = 1e-10)
      expect_lt(abs(fit$optimizer$gradient), 1e-4)
      expect_equal(unname(vcov(fit)['outcome.(Intercept)', 'outcome.(Intercept)']), .1, tolerance = 1e-6)
      expect_identical(fit$optimizer$control$parscale, scale)
    }
    if (backend == 'maxLik') {
      # maxLik's native Newton method uses code 1 for its gradient stopping rule.
      nr = panel_select(y ~ 1, NULL, data, 'id', 'poisson',
        start = c('outcome.(Intercept)' = log(2)), fixed = c(sd_re_outcome = 0),
        control = panel_control(optimizer = 'maxLik', method = 'NR',
          stationarity_warning = NULL), vcov = 'opg')
      expect_equal(nr$optimizer$code, 1)
      expect_true(nr$optimizer$converged)
      expect_equal(nr$loglik, maximum, tolerance = 1e-12)
    }
  })
}

test_that('the OPG quadratic is explicit, invariant to fixed scales and disableable', {
  score = matrix(c(1, 2, 1, 2, -1, 1), 3, 2,
    dimnames = list(NULL, c('beta', 'rho')))
  gradient = colSums(score); covariance = solve(crossprod(score))
  value = opg_stationarity(gradient, covariance)
  expect_identical(value$status, 'available')
  expect_equal(value$Q_OPG, 104/35, tolerance = 1e-14)
  expect_true(value$warning)
  expect_match(value$basis, '[Oo][Pp][Gg]')
  scale = c(beta = .1, rho = 10)
  transformed = opg_stationarity(gradient*scale, covariance/outer(scale, scale))
  expect_equal(transformed$Q_OPG, value$Q_OPG, tolerance = 1e-14)
  silent = opg_stationarity(gradient, covariance, tolerance = NULL)
  expect_equal(silent$Q_OPG, value$Q_OPG, tolerance = 1e-14)
  expect_false(silent$warning)
  invalid = covariance; invalid[1, 1] = -1
  expect_identical(opg_stationarity(gradient, invalid)$status, 'unavailable')
  expect_identical(opg_stationarity(numeric(), matrix(numeric(), 0, 0))$status, 'all_fixed')
  expect_error(panel_control(stationarity_warning = -1), 'stationarity_warning')
})

for (backend in c('optim', 'maxLik')) {
  test_that(paste('OPG warnings use supplied tolerances and historical diagnostics remain qualified:', backend), {
    if (backend == 'maxLik') skip_if_not_installed('maxLik')
    data = data.frame(id = 1:5, y = 0:4)
    formula = stats::as.formula('y ~ 1', env = baseenv())
    run = function(tolerance) panel_select(formula, NULL, data, 'id', 'poisson',
      start = c('outcome.(Intercept)' = log(.8)), fixed = c(sd_re_outcome = 0),
      control = panel_control(optimizer = backend, method = 'BFGS', parscale = 1,
        reltol = 2, stationarity_warning = tolerance), vcov = 'opg')
    # Deliberately loose objective stopping leaves a known nonzero score.
    expect_warning({ fit = run(1e-6) }, 'Q_OPG =')
    expect_equal(fit$optimizer$code, 0)
    expect_true(fit$optimizer$converged)
    mu = exp(unname(coef(fit)['outcome.(Intercept)']))
    expected = sum(data$y - mu)^2/sum((data$y - mu)^2)
    expect_equal(fit$optimizer$stationarity$Q_OPG, expected, tolerance = 1e-12)
    expect_true(fit$optimizer$stationarity$warning)
    expect_identical(fit$optimizer$stationarity$tolerance, 1e-6)
    for (tolerance in list(10, NULL)) {
      expect_warning({ quiet = run(tolerance) }, NA)
      expect_identical(quiet$optimizer$par, fit$optimizer$par)
      expect_identical(quiet$loglik, fit$loglik)
      expect_equal(quiet$optimizer$stationarity$Q_OPG, expected, tolerance = 1e-12)
      expect_false(quiet$optimizer$stationarity$warning)
      expect_identical(quiet$optimizer$stationarity$tolerance, tolerance)
    }
    stored = serialize(fit, NULL)
    value = check_stationarity(fit)
    expect_identical(value$status, 'unavailable')
    expect_match(value$reason, 'Observed information was not computed')
    expect_null(value$newton_step)
    expect_identical(value$opg$status, 'available')
    expect_equal(value$opg$Q_OPG, expected, tolerance = 1e-12)
    expect_identical(serialize(fit, NULL), stored)

    historical = fit
    historical$optimizer$stationarity = historical$optimizer$status = NULL
    historical$optimizer$rejections = historical$control$stationarity_warning = NULL
    value = check_stationarity(historical)
    expect_identical(value$status, 'unavailable')
    expect_identical(value$opg$status, 'available')
    expect_equal(value$opg$Q_OPG, expected, tolerance = 1e-12)
    expect_identical(value$opg$tolerance, .1)
    expect_null(value$optimizer_status)
    expect_null(value$rejections)
  })
}

test_that('all-fixed exact fits remain complete without an invented stopping test', {
  skip_if_not_installed('maxLik')
  data = data.frame(id = 1:3, y = c(1, 2, 0))
  fit = panel_select(y ~ 1, NULL, data, 'id', 'poisson',
    fixed = c('outcome.(Intercept)' = log(1.5), sd_re_outcome = 0),
    control = panel_control(optimizer = 'maxLik', maxit = 0), vcov = 'opg')
  expect_identical(fit$optimizer$status, 'all_fixed')
  expect_true(fit$optimizer$converged)
  expect_identical(fit$optimizer$stationarity$status, 'all_fixed')
  expect_false(fit$optimizer$stationarity$warning)
  expect_true(is.na(fit$optimizer$stationarity$Q_OPG))
  expect_equal(fit$loglik, sum(dpois(data$y, 1.5, log = TRUE)), tolerance = 1e-12)
})

test_that('accuracy-rejected physical split trials preserve the finite objective and score', {
  # The new default one-outcome identity resolves this singleton's tensor ridge.
  # Explicitly retain the original finite physical rule for rejection bookkeeping.
  physical = evaluate_split_binary
  local_mocked_bindings(evaluate_split_binary = function(theta, prepared, control,
      fixed = numeric(), score = TRUE) physical(theta, prepared, control, fixed,
        score, single_outcome = FALSE), .package = 'PanelSelect')
  data = data.frame(id = 1, d = 1, y = 1)
  prepared = prepare_data(y ~ 1, d ~ 1, data, 'id', 'probit', 'fail')
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0,
    sd_re_selection = 2.2, sd_re_outcome = 9, cor_re = 0, cor_error = 0)
  fixed = c(cor_error = 0)
  info = parameter_layout(prepared); info$fixed = info$name == 'cor_error'
  control = fit_control(panel_control(quadrature = 'gl', nodes = 20), vcov = 'opg')
  evaluator = fit_evaluator(theta, info, prepared, control, fixed)
  eta = parameter_to_optimizer(theta, info)
  before = evaluator$evaluate(eta, trial = FALSE)
  expect_equal(before$loglik, log(.25), tolerance = 1e-10)
  # The exact zero-index singleton derivative follows the arcsine probability.
  expected_rho = (2/pi)*2.2*9/sqrt((1 + 2.2^2)*(1 + 9^2))
  expect_equal(unname(before$optimizer_score[1, 'cor_re']), expected_rho, tolerance = 1e-9)

  trial = eta; trial['cor_re'] = atanh(.999)
  unavailable = evaluator$evaluate(trial)
  expect_false(unavailable$valid)
  expect_s3_class(unavailable$error, 'panelselect_split_accuracy_error')
  expect_identical(unavailable$loglik, -Inf)
  expect_identical(evaluator$rejected(), 1L)
  rejected = evaluator$rejections()
  expect_identical(rejected$total, 1L)
  expect_equal(sum(rejected$by_class), 1L)
  expect_equal(sum(rejected$by_reason), 1L)
  expect_true('panelselect_split_accuracy_error' %in% names(rejected$by_class))
  expect_true(any(grepl('likelihood|score|ambiguity', names(rejected$by_reason))))
  expect_true('panelselect_split_accuracy_error' %in% rejected$last$class)
  expect_identical(rejected$last$message, conditionMessage(unavailable$error))
  expect_identical(rejected$last$diagnostics, unavailable$error$diagnostics)
  expect_identical(evaluator$evaluate(trial), unavailable)
  expect_identical(evaluator$rejected(), 1L)
  expect_error(evaluator$evaluate(trial, trial = FALSE), class = 'panelselect_split_accuracy_error')
  after = evaluator$evaluate(eta, trial = FALSE)
  expect_identical(after$loglik, before$loglik)
  expect_identical(after$score, before$score)
  expect_identical(after$optimizer_score, before$optimizer_score)
  expect_identical(after$reference_optimizer_score, before$reference_optimizer_score)
  expect_identical(after$diagnostics$canonical_order, 20L)
  expect_identical(after$diagnostics$reference_order, 10L)
})
