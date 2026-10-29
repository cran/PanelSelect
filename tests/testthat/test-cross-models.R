cross_fixture = function(family, tau = .25) {
  panel = core_fixture(paste0("selection_", family))
  data = data.frame(id = c("b", "a", "b"), y = panel$prepared$y,
    d = panel$prepared$d, x = panel$prepared$x[, 2], w = panel$prepared$w[, 2])
  theta = panel$theta[!names(panel$theta) %in% c("sd_re_selection", "sd_re_outcome", "cor_re")]
  if ("cor_error" %in% names(theta)) theta["cor_error"] = tau
  prepared = prepare_data(y ~ x, d ~ w, data, NULL, family)
  list(data = data, theta = theta, prepared = prepared)
}

# Independent one-dimensional integrals; neither production CDF nor quadrature.
cross_reference = function(family, a, b, y, d, scale = .5, tau = .25) {
  integral = function(f, lower = -10, upper = 10)
    stats::integrate(f, lower, upper, rel.tol = 1e-10, abs.tol = 1e-12)$value
  if (!d) return(stats::pnorm(-a, log.p = TRUE))
  c = sqrt(1-tau^2)
  if (family == "probit") {
    q = 2*y-1
    return(log(integral(function(z) stats::dnorm(z) *
      stats::pnorm((q*b-q*tau*z)/c), upper = a)))
  }
  if (family == "gaussian") {
    z = (y-b)/scale
    return(stats::dnorm(z, log = TRUE)-log(scale) +
      stats::pnorm((a+tau*z)/c, log.p = TRUE))
  }
  if (family == "poisson")
    return(stats::pnorm(a, log.p = TRUE)+stats::dpois(y, exp(b), log = TRUE))
  log(integral(function(z) stats::dnorm(z)*stats::dpois(y, exp(b+scale*z))*
    stats::pnorm((a+tau*z)/c)))
}

test_that("cross-sectional row likelihoods and natural scores match independent mathematics", {
  control = panel_control(quadrature = "gh", observation_nodes = 41)
  for (family in c("probit", "gaussian", "pln")) for (tau in c(-.25, 0, .25)) {
    f = cross_fixture(family, tau)
    result = evaluate_likelihood(f$theta, f$prepared, control, numeric())
    a = drop(f$prepared$w %*% f$theta[startsWith(names(f$theta), "selection.")])
    b = drop(f$prepared$x %*% f$theta[startsWith(names(f$theta), "outcome.")])
    expected = vapply(seq_along(a), function(i) cross_reference(family, a[i], b[i],
      f$prepared$y[i], f$prepared$d[i], tau = tau), numeric(1))
    expect_equal(unname(result$loglik), expected, tolerance = 2e-9, info = family)
    info = parameter_layout(f$prepared)
    theta = f$theta[info$name]
    eta = parameter_to_optimizer(theta, info)
    transformed = colSums(result$score)[info$name]*parameter_jacobian(theta, info)
    for (step in c(1e-5, 3e-6)) {
      numerical = core_difference(function(theta) sum(evaluate_likelihood(theta,
        f$prepared, control, numeric(), score = FALSE)$loglik), f$theta, step)
      expect_equal(unname(colSums(result$score)), numerical, tolerance = 2e-6, info = family)
      numerical_eta = core_difference(function(value) sum(evaluate_likelihood(
        parameter_from_optimizer(value, theta, info), f$prepared, control,
        numeric(), score = FALSE)$loglik), eta, step)
      expect_equal(unname(transformed), numerical_eta, tolerance = 2e-6, info = family)
    }
    unused = f$prepared
    unused$y[unused$d == 0] = 999
    same = evaluate_likelihood(f$theta, unused, control, numeric())
    expect_equal(same$loglik, result$loglik)
    expect_equal(same$score, result$score)
  }
})

test_that("cross-sectional normalizations and fixed-zero count limits are exact", {
  a = -.3; b = .2; tau = .4; scale = .5
  binary = exp(cross_reference("probit", a, b, 0, 1, tau = tau)) +
    exp(cross_reference("probit", a, b, 1, 1, tau = tau))
  expect_equal(binary+pnorm(-a), 1, tolerance = 1e-11)
  normal = integrate(function(y) exp(vapply(y, function(value)
    cross_reference("gaussian", a, b, value, 1, scale, tau), numeric(1))), -8, 8)$value
  expect_equal(normal, pnorm(a), tolerance = 1e-8)
  pln = cross_fixture("pln")
  pln$theta[c("sd_error", "cor_error")] = 0
  control = panel_control(observation_nodes = 11)
  actual = evaluate_likelihood(pln$theta, pln$prepared, control,
    c(sd_error = 0, cor_error = 0))
  # The zero-scale PLN limit equals independent probit and Poisson regressions.
  prepared = pln$prepared
  iw = startsWith(names(pln$theta), "selection.")
  ix = startsWith(names(pln$theta), "outcome.")
  a = drop(prepared$w %*% pln$theta[iw])
  b = drop(prepared$x %*% pln$theta[ix])
  expected = vapply(seq_along(a), function(i) cross_reference("poisson", a[i], b[i],
    prepared$y[i], prepared$d[i]), numeric(1))
  sign = 2*prepared$d-1
  score_a = sign*dnorm(sign*a)/pnorm(sign*a)
  score_b = ifelse(prepared$d == 1, prepared$y-exp(b), 0)
  expected_score = cbind(prepared$w*score_a, prepared$x*score_b)
  expect_equal(unname(actual$loglik), expected, tolerance = 1e-13)
  expect_equal(unname(actual$score[, names(pln$theta)[iw | ix]]),
    unname(expected_score), tolerance = 1e-12)
})

test_that("cross predictions retain complete named natural Jacobians and covariance propagation", {
  for (family in c("probit", "gaussian", "pln")) {
    f = cross_fixture(family)
    fit = panel_select(y ~ x, d ~ w, f$data, id = NULL, family = family,
      fixed = f$theta, control = panel_control(observation_nodes = 21))
    V = diag(seq(.01, .02, length.out = length(f$theta)))
    dimnames(V) = list(names(coef(fit)), names(coef(fit)))
    targets = f$data[, c("x", "w")]
    for (type in c("response", "selection_probability", "joint", "selected_mean")) {
      result = predict(fit, targets, type = type, se.fit = TRUE, vcov = V)
      expect_identical(colnames(result$gradient), names(coef(fit)))
      numerical = vapply(seq_along(coef(fit)), function(j) {
        plus = minus = fit
        plus$coefficients[j] = plus$coefficients[j]+1e-5
        minus$coefficients[j] = minus$coefficients[j]-1e-5
        (predict(plus, targets, type = type)-predict(minus, targets, type = type))/2e-5
      }, numeric(nrow(targets)))
      expect_equal(unname(result$gradient), unname(numerical), tolerance = 2e-7)
      expect_equal(result$se.fit, sqrt(rowSums((result$gradient %*% V)*result$gradient)),
        ignore_attr = TRUE)
    }
    expect_error(predict(fit, targets, conditioning = "history", history = f$data),
      "cross|history", ignore.case = TRUE)
    expect_error(random_effects(fit), "cross|random", ignore.case = TRUE)
  }
})

test_that("cross model identity, criteria, serialization and updates use explicit structure", {
  f = cross_fixture("pln")
  f$theta[c("sd_error", "cor_error")] = 0
  fit = panel_select(y ~ x, d ~ w, f$data, id = NULL, family = "count",
    fixed = f$theta)
  expect_identical(fit$structure, "cross_sectional")
  expect_identical(fit$model_id, "probit_pln")
  expect_identical(model_structure(fit), "cross_sectional")
  expect_equal(attr(logLik(fit), "nobs"), nrow(f$data))
  expect_null(fit$sample)
  expect_identical(names(fit$loglik_by_unit), rownames(f$data))
  expect_equal(nobs(fit), nrow(f$data))
  count = fit
  count$parameter_info$fixed[1:2] = FALSE
  expect_equal(BIC(count), -2*fit$loglik+2*log(nrow(f$data)))
  expect_output(print(fit), "independent observations")
  expect_output(print(summary(fit)), "BIC \\(observations\\)")
  restored = unserialize(serialize(fit, NULL))
  expect_equal(coef(restored), coef(fit))
  expect_equal(predict(restored), predict(fit))
  call = update(restored, family = "pln", evaluate = FALSE)
  expect_true("id" %in% names(call))
  expect_null(call$id)
  expect_false(any(c("sd_re_selection", "sd_re_outcome", "cor_re") %in% names(call$start)))
  changed = update(restored, family = "pln")
  expect_identical(changed$model_id, "probit_pln")
  expect_equal(as.numeric(logLik(changed)), as.numeric(logLik(fit)), tolerance = 1e-12)
  coefficients = f$theta[!names(f$theta) %in% c("sd_error", "cor_error")]
  panel = update(restored, id = "id", family = "poisson", fixed = c(coefficients,
    sd_re_selection = 0, sd_re_outcome = 0, cor_re = 0))
  expect_identical(panel$structure, "panel")
  expect_identical(panel$model_id, "probit_re_poisson_re")
  expect_equal(as.numeric(logLik(panel)), as.numeric(logLik(fit)), tolerance = 1e-12)
  historical = panel
  historical$structure = historical$prepared$structure = NULL
  expect_identical(model_structure(historical), "panel")
  expect_equal(attr(logLik(historical), "nobs"), 2L)
  historical$prepared$id = NULL
  expect_error(update(historical, evaluate = FALSE), "missing metadata")
})
