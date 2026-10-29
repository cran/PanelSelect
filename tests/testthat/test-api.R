test_that("all seven models expose a consistent public fit", {
  for (model in core_models) {
    fixture = core_fixture(model)
    data = data.frame(id = c("b", "a", "b"), y = fixture$prepared$y,
      d = fixture$prepared$d, x = fixture$prepared$x[, 2], w = c(.5, -.3, .1))
    selection = if (fixture$prepared$has_selection) d ~ w else NULL
    # Independent value: panels "a" (row 2) and "b" (rows 1, 3) integrated by the
    # base-R reference with 40-node outer and 25-node inner rules. The Linear
    # route integrates its outcome effect exactly, so the reference needs 40 nodes.
    fine = panel_select(y ~ x, selection, data, "id", fixture$prepared$family,
      fixed = fixture$theta, control = panel_control(quadrature = "gh", nodes = 40L,
        observation_nodes = 25L))
    panel = function(rows) list(x = fixture$panel$x[rows, , drop = FALSE],
      w = fixture$panel$w[rows, , drop = FALSE], y = fixture$panel$y[rows],
      d = fixture$panel$d[rows])
    expected = ref_panel(model, fixture$reference, panel(2L), 40L, 25L)$loglik +
      ref_panel(model, fixture$reference, panel(c(1L, 3L)), 40L, 25L)$loglik
    expect_equal(as.numeric(logLik(fine)), expected, tolerance = 1e-9, info = model)
    for (method in c("gh", "agh")) {
      fit = panel_select(y ~ x,
        selection, data, "id", fixture$prepared$family, fixed = fixture$theta,
        control = panel_control(quadrature = method, nodes = 3, observation_nodes = 7))
      expect_s3_class(fit, "panelselect")
      expect_equal(coef(fit)[names(fixture$theta)], fixture$theta)
      expect_true(is.finite(as.numeric(logLik(fit))))
      expect_equal(attr(logLik(fit), "df"), 0L)
      expect_equal(nobs(fit), 3L)
      expect_equal(attr(logLik(fit), "nobs"), 2L)
      expect_true(all(vcov(fit) == 0))
      expect_true(all(is.na(summary(fit)$coefficients[, "se"])))
      expect_identical(fit$prepared$original_index, c(2L, 1L, 3L))
      expect_equal(model.matrix(fit)[, "x"], data$x, ignore_attr = TRUE)
      expect_equal(model.frame(fit)$y, data$y)
    }
  }
})

test_that("unused outcomes preserve selection information and missingness order", {
  data = data.frame(id = c(2, 1, 2, 1), d = c(1, 0, 1, 1),
    y = c(2, NA, 1, 0), x = c(.3, NA, -.2, .5), w = c(1, 2, 3, 4))
  prepared = prepare_data(y ~ x, d ~ w, data, "id", "poisson")
  expect_equal(prepared$n_rows, 4L)
  expect_equal(prepared$original_index, c(2L, 4L, 1L, 3L))
  expect_equal(prepared$x[prepared$d == 0, ], c("(Intercept)" = 0, x = 0))
  expect_equal(prepared$original_index[prepared$inverse_index], seq_len(4))
  data$w[3] = NA
  expect_error(prepare_data(y ~ x, d ~ w, data, "id", "poisson"), "row\\(s\\): 3")
  prepared = prepare_data(y ~ x, d ~ w, data, "id", "poisson", "exclude")
  expect_identical(prepared$keep, c(1L, 2L, 4L))
  expect_s3_class(prepared$na.action, "exclude")
  expect_equal(unname(as.integer(prepared$na.action)), 3L)
})

test_that("remembered designs retain factor levels in one-level history subsets", {
  data = data.frame(id = c(2, 1, 2, 1), y = c(0, 1, 2, 1),
                     category = c("a", "b", "a", "b"))
  prepared = prepare_data(y ~ factor(category), NULL, data, "id", "poisson")
  history = prepare_data(y ~ factor(category), NULL, data[1, , drop = FALSE],
    "id", "poisson", design = prepared$design)
  expect_identical(colnames(history$x), colnames(prepared$x))
  expect_equal(unname(history$x), matrix(c(1, 0), 1L))
  data$category[1] = "new"
  expect_error(prepare_data(y ~ factor(category), NULL, data[1, , drop = FALSE],
    "id", "poisson", design = prepared$design), "new level")
})

test_that("Poisson reduction reproduces GLM estimates and standard inference", {
  data = data.frame(id = rep(seq_len(6), each = 2),
    y = c(1, 2, 0, 3, 1, 0, 2, 1, 4, 0, 1, 3),
    x = c(-.8, .2, -.4, .7, .1, -.2, .5, -.3, .9, -.6, .3, .4),
    exposure = seq(.7, 1.3, length.out = 12))
  fit = panel_select(y ~ x + offset(log(exposure)), NULL, data, "id", "poisson",
                       fixed = c(sd_re_outcome = 0), vcov = "observed")
  reference = glm(y ~ x + offset(log(exposure)), data = data, family = poisson())
  expect_equal(unname(coef(fit)[1:2]), unname(coef(reference)), tolerance = 1e-5)
  expect_equal(unname(vcov(fit)[1:2, 1:2]), unname(vcov(reference)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)), tolerance = 1e-8)
  expect_equal(BIC(fit), -2*fit$loglik + 2*log(6))
  expect_equal(BIC(fit, unit = "observations"), -2*fit$loglik + 2*log(12))
  expect_equal(AIC(fit, fit)$AIC, rep(AIC(fit), 2))
  expect_equal(BIC(fit, fit, unit = "observations")$BIC,
               rep(BIC(fit, unit = "observations"), 2))
  expect_equal(BIC(fit, reference, unit = "observations")$BIC,
               rep(BIC(reference), 2), tolerance = 1e-8)
  expect_equal(unname(model.matrix(fit)), unname(model.matrix(reference)))
  expect_equal(formula(fit), y ~ x + offset(log(exposure)))
  expect_true(is.na(summary(fit)$coefficients["sd_re_outcome", "se"]))
  expect_equal(dim(confint(fit, parm = "outcome.x")), c(1L, 2L))
  call = update(fit, . ~ . - x, evaluate = FALSE)
  expect_equal(call$outcome, y ~ offset(log(exposure)))
  expect_error(vcov(fit, type = "automatic"), "arg")
})

for (backend in c("optim", "maxLik")) {
  test_that(paste("OPG-only fits compute explicitly requested observed inference:", backend), {
    if (backend == "maxLik") skip_if_not_installed("maxLik")
    data = data.frame(id = rep(seq_len(6), each = 2),
      y = c(1, 2, 0, 3, 1, 0, 2, 1, 4, 0, 1, 3),
      x = c(-.8, .2, -.4, .7, .1, -.2, .5, -.3, .9, -.6, .3, .4))
    fit = panel_select(y ~ x, NULL, data, "id", "poisson", fixed = c(sd_re_outcome = .3),
      control = panel_control(optimizer = backend, nodes = 3), vcov = "opg")
    expect_true(fit$optimizer$converged)
    expect_false(fit$optimizer$information_computed)
    expect_true(all(is.finite(vcov(fit, type = "observed"))))
    expect_true(all(is.finite(vcov(fit, type = "sandwich"))))
    if (backend == "maxLik") {
      reference = panel_select(y ~ x, NULL, data, "id", "poisson", fixed = c(sd_re_outcome = .3),
        control = panel_control(optimizer = "optim", nodes = 3), vcov = "opg")
      expect_equal(coef(reference), coef(fit), tolerance = 1e-5)
    }
  })
}

test_that("restrictions and optimizer scales remain distinct", {
  fixture = core_fixture("selection_gaussian")
  info = parameter_layout(fixture$prepared)
  fixed = resolve_fixed(c(sd_re_selection = 0), info, fixture$prepared)
  expect_equal(fixed[["cor_re"]], 0)
  expect_error(resolve_fixed(c(sd_re_selection = 0, cor_re = .2), info,
                              fixture$prepared), "cor_re")
  expect_error(resolve_fixed(c(sd_error = 0), info, fixture$prepared), "must be positive")
  expect_error(resolve_fixed(c(unknown = 0), info, fixture$prepared), "Unknown")
  info$fixed = !info$name %in% c("sd_re_outcome", "cor_re")
  theta = fixture$theta
  eta = parameter_to_optimizer(theta, info)
  expect_equal(unname(eta), c(log(.4), atanh(-.35)))
  expect_equal(parameter_from_optimizer(eta, theta, info), theta)
  information = matrix(c(3, 1, 1, 2), 2L)
  covariance = covariance_matrices(theta, info,
    matrix(0, 1L, length(theta)), information)$natural$observed
  jacobian = diag(c(.4, 1 - .35^2))
  expect_equal(unname(covariance[!info$fixed, !info$fixed]),
               jacobian %*% solve(information) %*% jacobian)
  expect_true(all(covariance[info$fixed, ] == 0))
})

test_that("obvious design and singleton probit identification failures are explicit", {
  data = data.frame(id = 1:4, y = c(0, 1, 0, 1), x = c(-.4, .3, .2, -.1))
  expect_error(panel_select(y ~ x, NULL, data, "id", "probit"), "one row per panel")
  expect_error(panel_select(y ~ x + I(2*x), NULL, data, "id", "probit",
                            fixed = c(sd_re_outcome = 0)), "not identified")
  fit = panel_select(y ~ x, NULL, data, "id", "probit", fixed = c(sd_re_outcome = 0), vcov = "observed")
  reference = glm(y ~ x, data = data, family = binomial("probit"))
  expect_equal(unname(coef(fit)[1:2]), unname(coef(reference)), tolerance = 1e-4)
})

test_that("observed information includes the second transform derivative", {
  fixture = core_fixture("poisson_re")
  info = parameter_layout(fixture$prepared)
  info$fixed = !info$name %in% c("outcome.(Intercept)", "sd_re_outcome")
  theta = fixture$theta
  fixed = theta[info$fixed]
  control = panel_control(nodes = 9)
  evaluator = fit_evaluator(theta, info, fixture$prepared, control, fixed)
  information = information_matrix(parameter_to_optimizer(theta, info), evaluator)
  free = which(!info$fixed)
  natural_score = function(value) colSums(evaluate_model(value, fixture$prepared,
    control, fixed)$score)[free]
  h = 1e-5
  natural_hessian = vapply(free, function(j) {
    plus = minus = theta
    plus[j] = plus[j] + h; minus[j] = minus[j] - h
    (natural_score(plus) - natural_score(minus))/(2*h)
  }, numeric(length(free)))
  jacobian = parameter_jacobian(theta, info)[free]
  second = c(0, theta[["sd_re_outcome"]])
  reference = -(natural_hessian*outer(jacobian, jacobian) +
                   diag(natural_score(theta)*second))
  expect_equal(unname(information), unname(reference), tolerance = 1e-5)
})

test_that("zero-column formulas have no phantom coefficients", {
  data = data.frame(id = c(2, 1, 2, 1), y = c(0, 1, 2, 1), d = c(0, 1, 1, 0))
  fit = panel_select(y ~ 0, NULL, data, "id", "poisson", fixed = c(sd_re_outcome = 0))
  expect_identical(names(coef(fit)), "sd_re_outcome")
  expect_equal(as.numeric(logLik(fit)), sum(dpois(data$y, 1, log = TRUE)))
  expect_equal(dim(model.matrix(fit)), c(4L, 0L))
  selected = panel_select(y ~ 0, d ~ 0, data, "id", "poisson",
    fixed = c(sd_re_selection = 0, sd_re_outcome = 0))
  expect_identical(names(coef(selected)), c("sd_re_selection", "sd_re_outcome", "cor_re"))
  expect_equal(as.numeric(logLik(selected)),
    nrow(data)*log(.5) + sum(dpois(data$y[data$d == 1], 1, log = TRUE)))
  expect_equal(dim(model.matrix(selected, equation = "selection")), c(4L, 0L))
})

test_that("complete starts preserve RNG and reject fixed-value conflicts", {
  data = data.frame(id = c(2, 1, 2, 1), y = c(0, 1, 2, 1), d = c(0, 1, 1, 0),
    category = factor(c(NA, "a", "b", NA)))
  prepared = prepare_data(y ~ category, d ~ 1, data, "id", "poisson")
  expect_true(all(is.finite(prepared$x)))
  info = parameter_layout(prepared)
  fixed = resolve_fixed(c(sd_re_outcome = 0), info, prepared)
  info$fixed = info$name %in% names(fixed)
  start = stats::setNames(rep(0, nrow(info)), info$name)
  start[info$type == "sd"] = 1
  start[names(fixed)] = fixed
  local_test_seed(22)
  seed = .Random.seed
  expect_equal(initial_parameters(prepared, info, start, fixed,
                                  fit_control(init = "norm")), start)
  expect_identical(.Random.seed, seed)
  expect_error(initial_parameters(prepared, info, c(cor_re = .2), fixed,
                                   panel_control()), "conflicts")
})
