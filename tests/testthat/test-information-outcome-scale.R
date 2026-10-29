# Gaussian outcome-unit invariance checked against independent references.
# The references below derive the selected density directly; they do not call
# package likelihoods, scores, parameter transforms or information-step helpers.

information_linear_data = function(offsets = FALSE) {
  local_test_seed(1301)
  n = 500L
  x = stats::rnorm(n); z = stats::rnorm(n); error = stats::rnorm(n)
  selection_error = .65*error + sqrt(1 - .65^2)*stats::rnorm(n)
  selection_offset = if (offsets) .15*sin(z) else rep(0, n)
  outcome_offset = if (offsets) .1*cos(x) else rep(0, n)
  d = as.integer(.2 + .4*x + .7*z + selection_offset + selection_error > 0)
  y = ifelse(d == 1, .3 + .5*x + outcome_offset + .7*error, NA_real_)
  data.frame(x = x, z = z, d = d, y = y,
    selection_offset = selection_offset, outcome_offset = outcome_offset)
}

information_linear_fit = function(data = information_linear_data(), start = NULL,
                                   fixed = numeric(), vcov = "opg", maxit = 200L) {
  outcome = stats::as.formula("y ~ x + offset(outcome_offset)", env = asNamespace("stats"))
  selection = stats::as.formula("d ~ x + z + offset(selection_offset)", env = asNamespace("stats"))
  panel_select(outcome, selection, data, NULL, "gaussian",
    start = start, fixed = fixed, init = "pooled", vcov = vcov,
    control = panel_control(quadrature = "gh", maxit = maxit, reltol = 1e-12))
}

# Full optimizer-coordinate Hessian and independent-unit scores. With v=log(sigma),
# t=atanh(tau), R=(y-b)/sigma and zeta=a*cosh(t)+R*sinh(t), the selected log density
# is log phi(R)-v+log Phi(zeta). Its log-Phi curvature is -M*(zeta+M).
information_linear_reference = function(fit) {
  data = fit$prepared; theta = coef(fit); parameters = names(theta)
  ia = match(paste0("selection.", colnames(data$w)), parameters)
  ib = match(paste0("outcome.", colnames(data$x)), parameters)
  iv = match("sd_error", parameters); it = match("cor_error", parameters)
  sigma = theta[["sd_error"]]; tau = theta[["cor_error"]]
  ch = 1/sqrt(1 - tau^2); sh = tau*ch
  a = drop(data$w %*% theta[ia]) + data$offset_selection
  b = drop(data$x %*% theta[ib]) + data$offset_outcome
  information = matrix(0, length(theta), length(theta), dimnames = list(parameters, parameters))
  score = matrix(0, nrow(data$x), length(theta), dimnames = list(NULL, parameters))
  for (i in seq_len(nrow(data$x))) {
    H = K = information*0; D = numeric(length(theta))
    if (data$d[i] == 0) {
      zeta = -a[i]
      M = exp(stats::dnorm(zeta, log = TRUE) - stats::pnorm(zeta, log.p = TRUE))
      H[ia, ia] = -M*(zeta + M)*tcrossprod(data$w[i, ])
      score[i, ia] = -M*data$w[i, ]
    } else {
      R = (data$y[i] - b[i])/sigma; zeta = a[i]*ch + R*sh
      M = exp(stats::dnorm(zeta, log = TRUE) - stats::pnorm(zeta, log.p = TRUE))
      D[ia] = data$w[i, ]*ch; D[ib] = -data$x[i, ]*sh/sigma
      D[iv] = -R*sh; D[it] = a[i]*sh + R*ch
      K[ia, it] = data$w[i, ]*sh
      K[ib, iv] = data$x[i, ]*sh/sigma; K[ib, it] = -data$x[i, ]*ch/sigma
      K[iv, it] = -R*ch
      K = K + t(K); K[iv, iv] = R*sh; K[it, it] = zeta
      H = -M*(zeta + M)*tcrossprod(D) + M*K
      H[ib, ib] = H[ib, ib] - tcrossprod(data$x[i, ])/sigma^2
      H[ib, iv] = H[ib, iv] - 2*R*data$x[i, ]/sigma
      H[iv, ib] = H[iv, ib] - 2*R*data$x[i, ]/sigma
      H[iv, iv] = H[iv, iv] - 2*R^2
      score[i, ] = M*D
      score[i, ib] = score[i, ib] + R*data$x[i, ]/sigma
      score[i, iv] = score[i, iv] + R^2 - 1
    }
    information = information - H
  }
  free = names(fit$optimizer$par)
  information = information[free, free, drop = FALSE]
  score = score[, free, drop = FALSE]
  observed = solve(information)
  list(information = information, score = score, observed = observed,
    sandwich = observed %*% crossprod(score) %*% observed)
}

# A separate scalar likelihood for Richardson differences, using only the
# observed selected/unselected density and the elementary coordinate transforms.
information_linear_loglik = function(eta, fit) {
  theta = coef(fit); theta[names(eta)] = eta
  sd = intersect(names(eta), "sd_error"); correlation = intersect(names(eta), "cor_error")
  theta[sd] = exp(theta[sd]); theta[correlation] = tanh(theta[correlation])
  data = fit$prepared; yes = data$d == 1
  a = drop(data$w %*% theta[paste0("selection.", colnames(data$w))]) + data$offset_selection
  b = drop(data$x %*% theta[paste0("outcome.", colnames(data$x))]) + data$offset_outcome
  sigma = theta[["sd_error"]]; tau = theta[["cor_error"]]
  R = (data$y[yes] - b[yes])/sigma
  sum(stats::pnorm(-a[!yes], log.p = TRUE)) + sum(stats::dnorm(R, log = TRUE) - log(sigma) +
    stats::pnorm((a[yes] + tau*R)/sqrt(1 - tau^2), log.p = TRUE))
}

information_linear_richardson = function(fit, h) {
  eta = fit$optimizer$par; p = length(eta)
  central = function(step) {
    H = matrix(0, p, p, dimnames = list(names(eta), names(eta)))
    value = function(point) information_linear_loglik(point, fit)
    middle = value(eta)
    for (j in seq_len(p)) for (k in seq_len(j)) {
      ej = ek = numeric(p); ej[j] = step; ek[k] = step
      H[j, k] = if (j == k) -(value(eta + ej) - 2*middle + value(eta - ej))/step^2 else
        -(value(eta + ej + ek) - value(eta + ej - ek) - value(eta - ej + ek) +
            value(eta - ej - ek))/(4*step^2)
      H[k, j] = H[j, k]
    }
    H
  }
  (4*central(h/2) - central(h))/3
}

expect_linear_information = function(fit, tolerance = 1e-5) {
  reference = information_linear_reference(fit)
  scale = outer(sqrt(diag(reference$information)), sqrt(diag(reference$information)))
  expect_lt(max(abs(fit$optimizer$information - reference$information)/scale), tolerance)
  for (type in c("observed", "sandwich")) {
    covariance = fit$optimizer$covariance[[type]]
    expect_true(all(is.finite(covariance)))
    expect_lt(max(abs(sqrt(diag(covariance)/diag(reference[[type]])) - 1)), tolerance)
  }
  invisible(reference)
}

test_that("Gaussian information has two independent references and respects outcome units", {
  base = information_linear_fit()
  reference = information_linear_reference(base)
  scale = outer(sqrt(diag(reference$information)), sqrt(diag(reference$information)))
  richardson = information_linear_richardson(base, .001)
  expect_lt(max(abs(richardson - reference$information)/scale), 1e-7)
  expect_lt(max(abs(richardson - information_linear_richardson(base, .002))/scale), 1e-7)
  expect_equal(information_linear_loglik(base$optimizer$par, base), base$loglik, tolerance = 1e-12)
  for (k in c(1, .001, .0001, .00001, 100)) {
    data = base$data; data$y = k*data$y
    theta = coef(base); outcome = startsWith(names(theta), "outcome.") | names(theta) == "sd_error"
    theta[outcome] = k*theta[outcome]
    expect_warning({ fit = information_linear_fit(data, theta, vcov = "observed", maxit = 0L) }, NA)
    expect_linear_information(fit)
    expect_equal(fit$loglik + sum(data$d)*log(k), base$loglik, tolerance = 1e-10)
  }
})

test_that("Gaussian information respects offsets, fixed parameters and joint regressor rescaling", {
  for (fixed in list(c(sd_error = .7, cor_error = .35), c(cor_error = 0))) {
    base = information_linear_fit(information_linear_data(TRUE), fixed = fixed)
    for (k in c(1e-4, 100)) {
      data = base$data; data$y = k*data$y; data$outcome_offset = k*data$outcome_offset
      data$x = 1000*data$x
      theta = coef(base); outcome = startsWith(names(theta), "outcome.") | names(theta) == "sd_error"
      theta[outcome] = k*theta[outcome]
      theta[c("outcome.x", "selection.x")] = theta[c("outcome.x", "selection.x")]/1000
      restrictions = theta[names(fixed)]
      expect_warning({ fit = information_linear_fit(data, theta, restrictions, "observed", 0L) }, NA)
      expect_linear_information(fit)
      expect_equal(fit$loglik + sum(data$d)*log(k), base$loglik, tolerance = 1e-10)
    }
  }
})

test_that("a normally optimized small-unit Gaussian fit uses corrected observed and sandwich information", {
  base = information_linear_fit()
  data = base$data; data$y = .001*data$y
  theta = coef(base); outcome = startsWith(names(theta), "outcome.") | names(theta) == "sd_error"
  theta[outcome] = .001*theta[outcome]
  theta["sd_error"] = 1.3*theta["sd_error"]
  expect_warning({ fit = information_linear_fit(data, theta, vcov = "observed") }, NA)
  expect_identical(fit$optimizer$status, "converged")
  expect_linear_information(fit)
})

test_that("information steps use the requested endpoint instead of the initial or cached SD", {
  data = information_linear_data()
  prepared = prepare_data(y ~ x, d ~ x + z, data, NULL, "gaussian")
  info = parameter_layout(prepared)
  theta = stats::setNames(c(.2, .4, .7, .3, .5, 2, .4), info$name)
  evaluator = fit_evaluator(theta, info, prepared, panel_control(quadrature = "gh"), numeric())
  initial = parameter_to_optimizer(theta, info)
  first = initial; first["sd_error"] = log(.002)
  second = initial; second["sd_error"] = log(.3)
  a = evaluator$information_steps(first); b = evaluator$information_steps(second)
  outcome = startsWith(names(a), "outcome.")
  expect_equal(unname(a["outcome.(Intercept)"]), 2e-6, tolerance = 1e-14)
  expect_equal(b[outcome]/a[outcome], stats::setNames(rep(150, sum(outcome)), names(a)[outcome]),
    tolerance = 1e-12)
  expect_identical(a[!outcome], b[!outcome])
  expect_identical(evaluator$information_steps(first), a)
  expect_identical(evaluator$evaluations(), 0L)
})

test_that("Gaussian on-demand covariance and stored information share the corrected calculation", {
  base = information_linear_fit()
  data = base$data; data$y = 1e-4*data$y
  theta = coef(base); scaled = startsWith(names(theta), "outcome.") | names(theta) == "sd_error"
  theta[scaled] = 1e-4*theta[scaled]
  fit = information_linear_fit(data, theta, maxit = 0L)
  saved = serialize(fit, NULL)
  expect_false(fit$optimizer$information_computed)
  reference = information_linear_reference(fit)
  jacobian = rep(1, length(theta)); names(jacobian) = names(theta)
  jacobian["sd_error"] = theta["sd_error"]
  jacobian["cor_error"] = 1 - theta["cor_error"]^2
  for (type in c("observed", "sandwich")) {
    expected = reference[[type]]*outer(jacobian, jacobian)
    expect_lt(max(abs(sqrt(diag(vcov(fit, type))/diag(expected)) - 1)), 1e-5)
    expect_lt(max(abs(summary(fit, vcov = type)$coefficients[, "se"]/sqrt(diag(expected)) - 1)), 1e-5)
  }
  stored = add_information(fit)
  expect_true(stored$optimizer$information_computed)
  expect_identical(stored$control$vcov, "opg")
  expect_identical(stored$covariance$opg, fit$covariance$opg)
  expect_identical(stored$coefficients, fit$coefficients)
  expect_linear_information(stored)
  expect_identical(add_information(stored), stored)
  expect_identical(serialize(fit, NULL), saved)
})

test_that("panel Gaussian information respects outcome units under GH, AGH and moving-cut GL", {
  data = data.frame(id = rep(1:6, each = 3), x = rep(c(-.6, .2, .7), 6),
    z = rep(c(.5, -.3, .1, .8, -.4, -.2), each = 3),
    d = rep(c(1, 0, 1, 1, 1, 0), 3), y = seq(-.5, .9, length.out = 18))
  data$y[data$d == 0] = NA_real_
  data$outcome_offset = .1*cos(data$x)
  data$selection_offset = .05*sin(data$z)
  theta = c("selection.(Intercept)" = .2, selection.z = .3,
    "outcome.(Intercept)" = .25, outcome.x = -.2, sd_re_selection = .35,
    sd_re_outcome = .4, sd_error = .7, cor_re = .3, cor_error = 0)
  fit_at = function(k, method) {
    current = data; current$y = k*current$y; current$outcome_offset = k*current$outcome_offset
    at = theta
    outcome = startsWith(names(at), "outcome.") | names(at) %in% c("sd_re_outcome", "sd_error")
    at[outcome] = k*at[outcome]
    fixed = at[!startsWith(names(at), "outcome.")]
    panel_select(y ~ x + offset(outcome_offset), d ~ z + offset(selection_offset),
      current, "id", "gaussian", start = at, fixed = fixed, vcov = "observed",
      control = panel_control(quadrature = method, nodes = 12L, maxit = 0L))
  }
  for (method in c("gh", "agh", "gl")) {
    expect_warning({ base = fit_at(1, method) }, NA)
    expect_warning({ scaled = fit_at(1e-4, method) }, NA)
    expected = base$optimizer$information/1e-8
    denominator = outer(sqrt(diag(expected)), sqrt(diag(expected)))
    expect_lt(max(abs(scaled$optimizer$information - expected)/denominator), 1e-7)
    expect_true(all(is.finite(scaled$optimizer$covariance$observed)))
    expect_lt(max(abs(sqrt(diag(scaled$optimizer$covariance$observed)/
      (1e-8*diag(base$optimizer$covariance$observed))) - 1)), 1e-7)
    if (method == "gl") {
      accuracy = attr(scaled$optimizer$information, "checked_accuracy")
      expect_identical(accuracy$status, "PASS")
      expect_equal(accuracy$native_optimizer_step,
        attr(base$optimizer$information, "checked_accuracy")$native_optimizer_step*1e-4,
        tolerance = 1e-12)
    }
  }
})

test_that("invalid information stencils retain an unavailable result without evaluating scores", {
  for (step in c(0, -1, NA_real_, Inf, .Machine$double.xmin)) {
    calls = 0L; reason = NULL
    evaluator = list(information_steps = function(eta) step,
      evaluate = function(...) { calls <<- calls + 1L; stop("unexpected score evaluation") })
    expect_warning({ H = information_matrix(c(beta = 1), evaluator,
      on_error = function(message) reason <<- message) }, "Observed information is unavailable")
    expect_true(all(is.na(H)))
    expect_match(reason, "positive and finite|not representable")
    expect_identical(calls, 0L)
  }
  evaluator = list(information_steps = function(eta) .Machine$double.eps,
    information_gate_tolerance = .01, evaluate = function(...) stop("unexpected score evaluation"))
  expect_warning({ H = information_matrix(c(beta = 1), evaluator) }, "not representable")
  expect_true(all(is.na(H)))
  evaluator$information_steps = function(eta) .Machine$double.xmax
  expect_warning({ H = information_matrix(c(beta = .Machine$double.xmax), evaluator) }, "not representable")
  expect_true(all(is.na(H)))
  expect_identical(information_matrix(numeric(), evaluator), matrix(numeric(), 0L, 0L))
})
