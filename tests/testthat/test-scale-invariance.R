# Scale invariance: rescaling a regressor or a Linear outcome rescales the fit,
# and rescaled optimizer coordinates leave the observed-information gate unchanged.

scale_data = function(seed = 3401, panels = 40L, waves = 4L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n), z = stats::rnorm(n))
  u = stats::rnorm(panels); v = .3*u + sqrt(1 - .09)*stats::rnorm(panels)
  data$d = as.integer(.3 + .4*data$x + .7*data$z + .6*u[data$id] + stats::rnorm(n) > 0)
  data$y = 1 + .5*data$x + .5*v[data$id] + stats::rnorm(n)
  data$y[data$d == 0] = NA
  data
}

# Retain the adaptive objective used to establish this optimizer-path comparison.
scale_fit = function(data) panel_select(y ~ x, d ~ x + z, data, "id", "linear",
  control = panel_control(quadrature = "agh", nodes = 6L), init = "pooled")

test_that("rescaled regressors and Linear outcomes give exactly rescaled fits", {
  data = scale_data()
  fit = scale_fit(data)
  small = data; small$x = small$x*1e-3; small$z = small$z*1e-3
  refit = scale_fit(small)
  slope = grepl("[.](x|z)$", names(coef(fit)))
  expected = coef(fit); expected[slope] = expected[slope]*1e3
  expect_equal(coef(refit), expected, tolerance = 1e-10)
  expect_equal(refit$loglik, fit$loglik, tolerance = 1e-12)
  expect_identical(refit$optimizer$counts, fit$optimizer$counts)
  outcome = grepl("^outcome[.]", names(coef(fit))) | names(coef(fit)) %in% c("sd_re_outcome", "sd_error")
  # 1e7 puts Linear outcome scores near 1e-7 of their regressor scale alone.
  for (scale in c(1e3, 1e7)) {
    large = data; large$y = large$y*scale
    refit = scale_fit(large)
    expected = coef(fit); expected[outcome] = expected[outcome]*scale
    expect_equal(coef(refit), expected, tolerance = 1e-10)
    # Each selected outcome density gains the Jacobian 1/scale.
    expect_equal(refit$loglik, fit$loglik - sum(data$d == 1)*log(scale), tolerance = 1e-12)
    expect_identical(refit$optimizer$counts, fit$optimizer$counts)
  }
})

test_that("coefficient scales follow their panel-score norms; scale coordinates stay capped", {
  # Small regressors give slope score norms below one, where a cap would bind.
  data = scale_data(); data$x = data$x*1e-3; data$z = data$z*1e-3
  fit = scale_fit(data)
  prepared = prepare_data(y ~ x, d ~ x + z, data, "id", "gaussian")
  info = parameter_layout(prepared)
  control = fit_control(panel_control(quadrature = "agh", nodes = 6L), init = "pooled")
  theta = initial_parameters(prepared, info, NULL, numeric(), control)
  score = fit_evaluator(theta, info, prepared, control, numeric())$evaluate(
    parameter_to_optimizer(theta, info), trial = FALSE)$optimizer_score
  norm = sqrt(colSums(score^2))
  coefficient = startsWith(colnames(score), "selection.") | startsWith(colnames(score), "outcome.")
  expected = pmax(sqrt(.Machine$double.eps), ifelse(coefficient, 1/norm, 1/pmax(1, norm)))
  expect_true(any(coefficient & norm < 1))
  expect_equal(unname(fit$optimizer$control$parscale), unname(expected), tolerance = 1e-14)
})

test_that("the observed-information gate is invariant to rescaled optimizer coordinates", {
  # Four nearby symmetric matrices (canonical, half step, predecessor and its half
  # step), then the same matrices for coordinates multiplied by 1/D: D H D.
  local_test_seed(3402)
  base = crossprod(matrix(stats::rnorm(16), 4L)) + diag(4)
  nearby = function(size) {
    change = matrix(stats::rnorm(16, sd = size), 4L)
    base + (change + t(change))/2
  }
  matrices = list(base, nearby(1e-5), nearby(3e-5), nearby(2e-5))
  scale = diag(10^c(-3, 0, 2, 5))
  rescaled = lapply(matrices, function(h) scale %*% h %*% scale)
  gate = function(m) information_accuracy_gate(m[[1L]], m[[2L]], m[[3L]], m[[4L]], .002)
  original = gate(matrices); scaled = gate(rescaled)
  # The statistic written out: entry changes over d_j d_k, d_j = sqrt(max(|a_jj|, |b_jj|)).
  standardized = function(a, b) max(abs(a - b)/outer(sqrt(pmax(abs(diag(a)), abs(diag(b)))),
    sqrt(pmax(abs(diag(a)), abs(diag(b))))))
  expect_equal(original$scaled_max_change, standardized(matrices[[1L]], matrices[[2L]]), tolerance = 1e-14)
  expect_equal(original$scaled_order_change, max(standardized(matrices[[1L]], matrices[[3L]]),
    standardized(matrices[[2L]], matrices[[4L]])), tolerance = 1e-14)
  expect_equal(scaled$scaled_max_change, original$scaled_max_change, tolerance = 1e-10)
  expect_equal(scaled$scaled_order_change, original$scaled_order_change, tolerance = 1e-10)
  expect_identical(scaled$status, original$status)
  expect_identical(original$status, "PASS")
})

test_that("integration score floors carry regressor units", {
  data = scale_data()
  prepared = prepare_data(y ~ x, d ~ x + z, data, "id", "gaussian")
  names = parameter_layout(prepared)$name
  theta = stats::setNames(numeric(length(names)), names)
  floors = score_floors(prepared, theta)
  rms = function(v) { v = v[is.finite(v) & v != 0]; sqrt(mean(v^2)) }
  expect_equal(floors[["selection.z"]], rms(data$z))
  expect_equal(floors[["outcome.x"]], rms(data$x[data$d == 1]))
  scale = !grepl("^(selection|outcome)[.]", names)
  expect_identical(unname(floors[c("selection.(Intercept)", "outcome.(Intercept)", names[scale])]),
    rep(1, 2L + sum(scale)))
  # An evaluator workspace keeps the first result for its fixed data.
  workspace = new.env(parent = emptyenv())
  first = score_floors(prepared, theta, workspace)
  expect_identical(score_floors(NULL, theta, workspace), first)
})

test_that("staged and split integration checks do not depend on covariate units", {
  simulate = function(family, scale) {
    local_test_seed(5301)
    panels = 12L; waves = 3L; id = rep(seq_len(panels), each = waves); n = panels*waves
    x = stats::rnorm(n); w = stats::rnorm(n)
    u = stats::rnorm(panels)[id]; v = .4*u + sqrt(.84)*stats::rnorm(panels)[id]
    d = as.integer(.3 + .5*w + .5*u + stats::rnorm(n) > 0)
    index = .2 + .4*x + .6*v
    y = if (family == "poisson") stats::rpois(n, exp(index)) else as.integer(index + stats::rnorm(n) > 0)
    y[d == 0] = NA
    data.frame(id = id, x = x*scale, w = w*scale, d = d, y = y)
  }
  evaluate = function(family, scale, fixed) {
    unscaled = prepare_data(y ~ x, d ~ x + w, simulate(family, 1), "id", family, "fail")
    theta = initial_parameters(unscaled, parameter_layout(unscaled), NULL, fixed,
      fit_control(init = "pooled"))
    slope = grepl("[.](x|w)$", names(theta)); theta[slope] = theta[slope]/scale
    prepared = prepare_data(y ~ x, d ~ x + w, simulate(family, scale), "id", family, "fail")
    # Selection Poisson is staged under GL; probit with cor_error fixed at 0 takes moving cuts.
    control = staged_test_budget(panel_control(quadrature = "gl"))
    staged = uses_staged_integration(prepared, control, fixed)
    if (staged) control$integration_mesh = prepare_staged_integration(theta, prepared, control, fixed)
    value = evaluate_model(theta, prepared, control, fixed)
    list(loglik = value$loglik, score = value$score, slope = slope,
      nodes = if (staged) vapply(control$integration_mesh$outer, staged_node_count, numeric(1)))
  }
  for (case in list(list("poisson", numeric()), list("probit", c(cor_error = 0)))) {
    base = evaluate(case[[1]], 1, case[[2]])
    # With absolute score budgets the staged check required a new mesh at 1e4.
    scaled = evaluate(case[[1]], 1e4, case[[2]])
    expect_identical(scaled$nodes, base$nodes)
    expect_equal(scaled$loglik, base$loglik, tolerance = 1e-12)
    expected = base$score; expected[, base$slope] = expected[, base$slope]*1e4
    expect_equal(scaled$score, expected, tolerance = 1e-9)
  }
})
