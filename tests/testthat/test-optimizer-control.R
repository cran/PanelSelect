# Behavioral coverage ported from master c9f2642 to the shared version 2 engine.
control_fixture = function(binary) {
  data = data.frame(id = rep(1:8, each = 3), x = rep(c(-.8, .1, .7), 8),
    d = rep(c(1, 0, 1, 1, 1, 0), 4), y = rep(c(0, NA, 1, 1, 0, NA), 4))
  if (!binary) data$y = data$y + seq_len(nrow(data)) / 20
  data
}

# The optimizer options of a binary selection fit from zero starts with both
# correlations free, captured where optimization would begin.
capture_options = function(control) {
  observed = NULL
  local_mocked_bindings(optimize_model = function(eta, evaluator, control, options) {
    observed <<- options
    stop("optimizer captured")
  })
  expect_error(panel_select(y ~ x, d ~ x, control_fixture(TRUE), "id", "probit", init = "zero",
    control = control), "optimizer captured")
  observed
}

test_that("the optimizer minimizes the negative log likelihood without fnscale", {
  evaluator = list(evaluate = function(eta)
    list(loglik = -eta^2, optimizer_score = matrix(-2*eta, 1)))
  called = FALSE
  local_mocked_bindings(optim = function(par, fn, gr, method, control) {
    called <<- TRUE
    expect_equal(fn(par), 9)
    expect_equal(gr(par), 6)
    expect_null(control$fnscale)
    expect_equal(control$parscale, 2)
    list(par = 0, convergence = 0)
  }, .package = "stats")
  expect_equal(optimize_model(3, evaluator, panel_control(parscale = 2))$par, 0)
  expect_true(called)
})

test_that('automatic binary BFGS tuning is scoped and respects explicit controls', {
  prepared = list(model_id = 'probit_re_probit_re', n_panels = 400L)
  info = data.frame(name = c('outcome.x', 'cor_error'), fixed = c(TRUE, FALSE))
  control = panel_control(quadrature = 'gh')
  options = optimizer_options(validate_panel_control(control), prepared, info)
  expect_identical(options$parscale, 1/20)
  expect_identical(options$reltol, 1e-12)
  expect_identical(control$reltol, sqrt(.Machine$double.eps))
  for (explicit in list(panel_control(quadrature = 'gh', reltol = sqrt(.Machine$double.eps)),
      panel_control(quadrature = 'gh', reltol = 1e-10),
      list(quadrature = 'gh', reltol = sqrt(.Machine$double.eps)),
      panel_control(quadrature = 'gh', parscale = 2))) {
    actual = optimizer_options(validate_panel_control(explicit), prepared, info)
    expected = optimizer_options(validate_panel_control(explicit))
    expect_identical(actual, expected)
  }
  changed = control; changed$reltol = 1e-8
  expect_null(optimizer_options(validate_panel_control(changed), prepared, info)$parscale)
  for (field in c('quadrature', 'optimizer', 'method')) {
    other = control
    other[[field]] = switch(field, quadrature = 'agh', optimizer = 'maxLik', method = 'CG')
    expect_null(optimizer_options(other, prepared, info)$parscale)
  }
  for (model in c('probit_re', 'probit_probit', 'probit_re_linear_re')) {
    other = prepared; other$model_id = model
    expect_null(optimizer_options(control, other, info)$parscale)
  }
  info$fixed[] = TRUE
  expect_null(optimizer_options(control, prepared, info)$parscale)
  for (extra in list(list(), list(maxit = 17))) {
    auto = capture_options(do.call(panel_control, c(list(quadrature = 'gh', nodes = 3L), extra)))
    expect_equal(auto$parscale, rep(1/sqrt(8), 8L))
    expect_identical(auto$reltol, 1e-12)
  }
  for (extra in list(list(reltol = sqrt(.Machine$double.eps)), list(parscale = 2))) {
    explicit = capture_options(do.call(panel_control, c(list(quadrature = 'gh', nodes = 3L), extra)))
    if (is.null(extra$parscale)) {
      # An explicit tolerance disables the older binary tuning pair, but
      # retains the new independent initial-score scaling.
      prepared = prepare_data(y ~ x, d ~ x, control_fixture(TRUE), 'id', 'probit', 'fail')
      layout = parameter_layout(prepared)
      layout$fixed = FALSE
      theta = initial_parameters(prepared, layout, NULL, numeric(), fit_control(init = 'zero'))
      value = evaluate_model(theta, prepared, panel_control(quadrature = 'gh', nodes = 3L), numeric())
      score = sweep(value$score, 2L, parameter_jacobian(theta, layout), '*')
      # Coefficient scales follow their panel-score norms unless the
      # score is structurally zero (here outcome.(Intercept) at the zero start);
      # log-SD and atanh-correlation scales keep the cap at one.
      norm = sqrt(colSums(score^2))
      reference = design_scales(prepared)[colnames(score)]*sqrt(prepared$n_panels)
      uncapped = !is.na(reference) & norm > 1e-6*reference
      expect_equal(explicit$parscale, unname(ifelse(uncapped, 1/norm, 1/pmax(1, norm))))
    } else expect_equal(explicit$parscale, rep(extra$parscale, 8L))
    expect_equal(explicit$reltol, extra$reltol %||% sqrt(.Machine$double.eps))
  }
})

test_that('binary tuning is recomputed after serialization, update and refinement', {
  data = control_fixture(TRUE)
  fixed = c('selection.(Intercept)' = .1, selection.x = -.1,
    'outcome.(Intercept)' = .2, outcome.x = .1, sd_re_selection = .4, sd_re_outcome = .3, cor_re = .2)
  control = panel_control(quadrature = 'gh', nodes = 3L, maxit = 0L)
  fit = panel_select(y ~ x, d ~ x, data, 'id', 'probit',
    start = c(cor_error = .1), fixed = fixed, control = control)
  saved = unserialize(serialize(fit, NULL))
  doubled = rbind(data, transform(data, id = id+8))
  updated = update(saved, data = doubled)
  refined = check_quadrature(saved, nodes = 4L, refit = TRUE)$refit
  for (object in list(saved, updated, refined)) {
    expect_identical(attr(object$control, 'reltol_supplied'), FALSE)
    expect_identical(object$optimizer$control$reltol, 1e-12)
    expect_equal(object$optimizer$control$parscale, 1/sqrt(object$prepared$n_panels))
  }
  explicit = update(saved, control = panel_control(quadrature = 'gh',
    nodes = 3L, maxit = 0L, reltol = sqrt(.Machine$double.eps)))
  expect_null(explicit$control$parscale)
  score = explicit$score[, 'cor_error'] * (1 - explicit$coefficients[['cor_error']]^2)
  expect_equal(explicit$optimizer$control$parscale,
    unname(1/max(1, sqrt(sum(score^2)))))
  expect_true(attr(explicit$control, 'reltol_supplied'))
})
