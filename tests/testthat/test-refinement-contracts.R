# Regression contracts for the readability/inference refinement. The Poisson
# score and information below come directly from the scalar log likelihood.
refinement_poisson_fit = function(vcov = "opg", ...) {
  data = data.frame(id = 1:5, y = 0:4, x = seq(-1, 1, length.out = 5))
  panel_select(y ~ 1, NULL, data, "id", "poisson",
    start = c("outcome.(Intercept)" = log(1.5)), fixed = c(sd_re_outcome = 0),
    control = panel_control(maxit = 0, ...), vcov = vcov)
}

test_that("Q diagnostics distinguish OPG and observed curvature without square roots", {
  fit = refinement_poisson_fit("observed")
  saved = serialize(fit, NULL)
  actual = check_stationarity(fit)
  # g = 10 - 5*1.5; B = sum((0:4 - 1.5)^2); A = 5*1.5.
  expect_equal(actual$Q_OPG, 5/9, tolerance = 1e-12)
  expect_equal(actual$Q_Hessian, 5/6, tolerance = 1e-6)
  expect_equal(actual$Q_Hessian, sum(actual$score*actual$newton_step), tolerance = 1e-12)
  expect_equal(fit$optimizer$stationarity$Q_OPG, actual$Q_OPG)
  expect_equal(fit$optimizer$stationarity$Q_Hessian, actual$Q_Hessian)
  expect_false(any(c("score_norm", "quadratic") %in% names(actual$opg)))
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "Q_OPG =")
  expect_match(paste(capture.output(summary(fit)), collapse = "\n"), "Q_Hessian =")
  expect_identical(serialize(fit, NULL), saved)
})

test_that("the Q_OPG warning uses a strict quadratic threshold", {
  covariance = matrix(1, 1, 1, dimnames = list("b", "b"))
  value = opg_stationarity(c(b = 1), covariance, tolerance = 1)
  expect_equal(value$Q_OPG, 1)
  expect_false(value$warning)
  expect_true(opg_stationarity(c(b = 1), covariance, tolerance = 1 - 1e-8)$warning)
  expect_false(opg_stationarity(c(b = 1), covariance, tolerance = 1 + 1e-8)$warning)
  # This lies between the old default Q cutoff .01 and the new cutoff .1.
  value = opg_stationarity(c(b = .2), covariance)
  expect_equal(value$Q_OPG, .04)
  expect_false(value$warning)
  expect_true(opg_stationarity(c(b = .4), covariance)$warning)
  expect_false(opg_stationarity(c(b = 2), covariance, tolerance = NULL)$warning)
  expect_identical(panel_control()$stationarity_warning, .1)
  tenth = covariance/10
  expect_equal(opg_stationarity(c(b = 1), tenth)$Q_OPG, .1)
  expect_false(opg_stationarity(c(b = 1), tenth)$warning)
  expect_true(opg_stationarity(c(b = 1 + 1e-8), tenth)$warning)
  expect_false(opg_stationarity(c(b = 1 - 1e-8), tenth)$warning)
  overflow = opg_stationarity(c(b = 1e200), covariance)
  expect_identical(overflow$status, "unavailable")
  expect_true(is.na(overflow$Q_OPG))
})

test_that("Hessian Q is unavailable for invalid stored information or inverse", {
  fit = refinement_poisson_fit("observed")
  for (problem in c("missing", "singular", "indefinite", "nonfinite", "inverse")) {
    bad = fit
    if (problem == "missing") bad$optimizer$information_computed = FALSE
    if (problem == "singular") bad$optimizer$information[1, 1] = 0
    if (problem == "indefinite") bad$optimizer$information[1, 1] = -1
    if (problem == "nonfinite") bad$optimizer$information[1, 1] = Inf
    if (problem == "inverse") bad$optimizer$covariance$observed[1, 1] = -1
    checked = check_stationarity(bad)
    expect_true(is.na(checked$Q_Hessian), info = problem)
    expect_identical(checked$status, "unavailable", info = problem)
    expect_true(nzchar(checked$reason), info = problem)
    expect_equal(checked$Q_OPG, 5/9)
  }
  all_fixed = panel_select(y ~ 1, NULL, data.frame(id = 1:3, y = 0:2), "id", "poisson",
    fixed = c("outcome.(Intercept)" = 0, sd_re_outcome = 0))
  checked = check_stationarity(all_fixed)
  expect_identical(checked$status, "all_fixed")
  expect_true(is.na(checked$Q_OPG))
  expect_true(is.na(checked$Q_Hessian))
})

test_that("summaries reuse requested information and leave the supplied fit unchanged", {
  fit = refinement_poisson_fit()
  saved = serialize(fit, NULL)
  calls = 0L
  original = information_matrix
  local_mocked_bindings(information_matrix = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "PanelSelect")
  expect_true(is.na(check_stationarity(fit)$Q_Hessian))
  invisible(capture.output(print(fit)))
  expect_true(is.na(summary(fit)$stationarity$Q_Hessian))
  expect_identical(calls, 0L)
  for (type in c("observed", "sandwich")) {
    before = calls
    answer = summary(fit, vcov = type)
    expect_identical(calls, before + 1L)
    expect_equal(answer$stationarity$Q_Hessian, 5/6, tolerance = 1e-6)
    expect_identical(serialize(fit, NULL), saved)
  }
  enriched = add_information(fit)
  before = calls
  expect_equal(check_stationarity(enriched)$Q_Hessian, 5/6, tolerance = 1e-6)
  invisible(summary(enriched, vcov = "sandwich"))
  expect_identical(calls, before)
  expect_identical(add_information(enriched), enriched)
})

test_that("saved native fit displays recompute the new Q flag from stored quantities", {
  fit = refinement_poisson_fit()
  fit$optimizer$gradient[] = .2
  fit$optimizer$covariance$opg[,] = 1
  fit$control$quadrature = "agh"
  fit$optimizer$stationarity = list(status = "available", score_norm = .2,
    quadratic = .04, tolerance = .1, warning = TRUE, reason = NULL, basis = "panel_opg")
  saved = serialize(fit, NULL)
  checked = check_stationarity(fit)
  expect_equal(checked$Q_OPG, .04)
  expect_false(checked$opg$warning)
  expect_false(summary(fit)$stationarity$warning)
  expect_false(grepl("substantial|OPG-scaled score", paste(capture.output(print(fit)), collapse = "\n")))
  expect_identical(serialize(fit, NULL), saved)
  expect_identical(fit$control$quadrature, "agh")
  fit$control$stationarity_warning = .02
  expect_true(check_stationarity(fit)$opg$warning)
  fit$control$stationarity_warning = NULL
  expect_false(check_stationarity(fit)$opg$warning)
  old_summary = summary(fit)
  old_summary$stationarity = list(quadratic = Inf, score_norm = Inf,
    tolerance = .1, warning = TRUE)
  printed = NULL
  invisible(capture.output({ printed = print(old_summary) }))
  expect_true(is.na(printed$stationarity$Q_OPG))
  expect_false(printed$stationarity$warning)
})

test_that("update validates supplied starts before warm-start merging", {
  fit = refinement_poisson_fit()
  bad = list(42, matrix(.2, 1, dimnames = list(NULL, "outcome.(Intercept)")),
    c("outcome.(Intercept)" = .2, "outcome.(Intercept)" = .4),
    setNames(.2, ""), setNames(.2, NA_character_), c("outcome.(Intercept)" = Inf),
    c("outcome.(Intercept)" = NA_real_), c(unknown = .1), c(sd_re_outcome = -1))
  for (start in bad) for (evaluate in c(FALSE, TRUE))
    expect_error(update(fit, start = start, evaluate = evaluate), "start|parameter|name|finite|SD")
  partial = update(fit, start = c("outcome.(Intercept)" = .7), evaluate = FALSE)
  expect_equal(partial$start[["outcome.(Intercept)"]], .7)
  expect_equal(coef(update(fit, start = c("outcome.(Intercept)" = .7)))[["outcome.(Intercept)"]], .7)
  expect_null(update(fit, start = NULL, evaluate = FALSE)$start)
  expect_equal(update(fit, start = numeric(), evaluate = FALSE)$start, coef(fit))
  changed = update(fit, . ~ . + x, start = c(outcome.x = .3), evaluate = FALSE)
  expect_equal(changed$start[["outcome.x"]], .3)
  freed = update(fit, fixed = numeric(), evaluate = FALSE)
  expect_false("sd_re_outcome" %in% names(freed$start))
  expect_error(update(update(fit, . ~ . + x), . ~ . - x,
    start = c(outcome.x = .2), evaluate = FALSE), "Unknown|unknown")
  data = data.frame(y = c(0, 1, 0, 1, NA, NA), d = c(1, 1, 1, 1, 0, 0))
  binary = panel_select(y ~ 1, d ~ 1, data, id = NULL, family = "probit",
    fixed = c("outcome.(Intercept)" = 0), control = panel_control(maxit = 0))
  for (correlation in c(-1, 1, 1.01)) for (evaluate in c(FALSE, TRUE))
    expect_error(update(binary, start = c(cor_error = correlation), evaluate = evaluate),
      "correlation|Correlation|strictly")
})

test_that("default GH orders and stored rules have distinct contracts", {
  control = panel_control()
  expect_identical(control$quadrature, "gh")
  expect_equal(control$nodes, 10)
  expect_identical(control$observation_quadrature, "gh")
  expect_equal(control$observation_nodes, 20)
  data = data.frame(id = rep(1:5, each = 2), y = rep(0:4, each = 2))
  fit_with_rule = function(control) panel_select(y ~ 1, NULL, data, "id", "poisson",
    start = c("outcome.(Intercept)" = .3, sd_re_outcome = .6), control = control)
  default = fit_with_rule(panel_control(maxit = 0))
  explicit = fit_with_rule(panel_control(quadrature = "gh", nodes = 10, maxit = 0))
  expect_identical(coef(default), coef(explicit))
  expect_identical(logLik(default), logLik(explicit))
  for (method in c("agh", "gl")) {
    old = refinement_poisson_fit(quadrature = method)
    call = update(old, evaluate = FALSE)
    expect_identical(call$control$quadrature, method)
    expect_identical(update(old)$control$quadrature, method)
  }
})

test_that("BIC sample size and GL history restrictions remain unchanged", {
  data = data.frame(id = rep(1:5, each = 2), y = rep(0:4, each = 2))
  fit = panel_select(y ~ 1, NULL, data, "id", "poisson", fixed = c(sd_re_outcome = 0))
  k = attr(logLik(fit), "df")
  expect_equal(nobs(fit), 10)
  expect_equal(attr(logLik(fit), "nobs"), 5)
  expect_equal(BIC(fit), -2*as.numeric(logLik(fit)) + log(5)*k)
  expect_equal(BIC(fit, unit = "observations"), -2*as.numeric(logLik(fit)) + log(10)*k)
  fit$control$quadrature = "gl"
  expect_error(predict(fit, conditioning = "history", history = data), "unavailable for GL")
  expect_error(random_effects(fit, data), "GL")
})
