# Fitting-interface inputs: the outcome argument, the .outcome shorthand,
# rejected dots, the Linear family aliases, the explicit cross-sectional id,
# subset rows, specification messages and control defaults.

# Twenty panels of four waves with interior estimates (no boundary correlation).
interface_data = function(seed = 1101, panels = 20L, waves = 4L, family = 'poisson') {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + stats::rnorm(n) > 0)
  index = .2 + .4*data$x + .6*v[data$id]
  data$y = if (family == 'poisson') stats::rpois(n, exp(index)) else index + .5*stats::rnorm(n)
  data$y[data$d == 0] = NA
  data
}

# Tight stopping keeps these small fits clear of the OPG-score screen; the fits
# below start from pooled regressions (init = 'pooled').
interface_control = function(...) panel_control(nodes = 4L, reltol = 1e-12, ...)

expect_same_fit = function(actual, expected) {
  expect_identical(actual$coefficients, expected$coefficients)
  expect_identical(actual$loglik, expected$loglik)
  expect_identical(actual$score, expected$score)
}

test_that('outcome is the first argument and the formula alias is gone', {
  data = interface_data()
  control = interface_control()
  positional = panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', control = control,
    init = 'pooled')
  named = panel_select(outcome = y ~ x, selection = d ~ x + w, data = data, id = 'id',
    family = 'poisson', control = control, init = 'pooled')
  expect_same_fit(named, positional)
  expect_identical(names(positional$call)[2L], 'outcome')
  expect_equal(formula(positional), y ~ x)
  expect_error(panel_select(formula = y ~ x, selection = d ~ x + w, data = data, id = 'id',
    family = 'poisson'), 'unused argument')
})

test_that('dots are rejected in both equations and in standalone models', {
  data = interface_data()
  expect_error(panel_select(y ~ ., d ~ x + w, data, 'id', 'poisson'), 'not supported in the outcome')
  expect_error(panel_select(y ~ x, d ~ . - y, data, 'id', 'poisson'), 'not supported in the selection')
  observed = data[data$d == 1, ]
  expect_error(panel_select(y ~ ., NULL, observed, 'id', 'poisson'), 'not supported in the outcome')
  # A stored design is reused as fitted: an old object's formula is not re-parsed.
  fit = panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', control = interface_control(),
    init = 'pooled')
  fit$prepared$design$outcome$formula = y ~ .
  expect_equal(predict(fit, data[1:3, ], conditioning = 'history', history = data),
    predict(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', control = interface_control(),
      init = 'pooled'), data[1:3, ], conditioning = 'history', history = data))
})

test_that('.outcome stands for the outcome regressors in the selection equation', {
  data = interface_data()
  data$f = factor(rep(c('a', 'b', 'c'), length.out = nrow(data)))
  data$`my var` = seq(-1, 1, length.out = nrow(data))
  data$o = seq(.1, .3, length.out = nrow(data))
  control = interface_control()
  fit = function(outcome, selection) panel_select(outcome, selection, data, 'id', 'poisson',
    control = control, init = 'pooled')
  short = fit(y ~ x, d ~ .outcome + w)
  expect_same_fit(short, fit(y ~ x, d ~ x + w))
  expect_equal(formula(short, 'selection'), d ~ (x) + w, ignore_attr = TRUE)
  # Richer equations compare their prepared selection designs and offsets.
  same_design = function(outcome, shorthand, explicit) {
    a = prepare_data(outcome, shorthand, data, 'id', 'poisson')
    b = prepare_data(outcome, explicit, data, 'id', 'poisson')
    expect_identical(a$w, b$w)
    expect_identical(a$offset_selection, b$offset_selection)
    expect_identical(a$design$selection$shorthand, shorthand)
  }
  same_design(y ~ poly(x, 2) + f*`my var`, d ~ .outcome + w,
    d ~ poly(x, 2) + f + `my var` + f:`my var` + w)
  # Offsets and the outcome's intercept setting are not carried over.
  same_design(y ~ 0 + x + offset(o), d ~ .outcome + w, d ~ x + w)
  same_design(y ~ x, d ~ w + .outcome:w, d ~ w + x:w)
  expect_error(fit(y ~ x, d ~ log(.outcome) + w), 'not an argument of log')
  expect_error(fit(y ~ x, .outcome ~ w), 'cannot be the selection response')
  expect_error(fit(y ~ .outcome, d ~ w), 'only in the selection formula')
  expect_error(fit(y ~ 1, d ~ .outcome + w), 'outcome formula has none')
})

test_that('linear aliases give the same fit and id must be written out', {
  data = interface_data(family = 'gaussian')
  # Fits at fixed starts (maxit = 0) compare the likelihood without optimizing.
  control = interface_control(maxit = 0L)
  theta = c('selection.(Intercept)' = .3, selection.w = .6, 'outcome.(Intercept)' = .2,
    outcome.x = .4, sd_re_selection = .8, sd_re_outcome = .6, sd_error = .5, cor_re = .4, cor_error = .2)
  gaussian = panel_select(y ~ x, d ~ w, data, 'id', 'gaussian', start = theta, control = control,
    init = 'pooled')
  linear = panel_select(y ~ x, d ~ w, data, 'id', 'linear', start = theta, control = control,
    init = 'pooled')
  expect_same_fit(linear, gaussian)
  expect_identical(linear$prepared$family, 'gaussian')
  expect_same_fit(update(gaussian, family = 'linear'), update(gaussian))
  continuous = panel_select(y ~ x, d ~ w, data, 'id', 'continuous', start = theta,
    control = control, init = 'pooled')
  expect_same_fit(continuous, linear)
  expect_same_fit(update(gaussian, family = 'continuous'), update(gaussian))
  expect_error(panel_select(y ~ x, d ~ w, data, family = 'gaussian'),
    'id = NULL for independent cross-sectional observations')
  cross = data; cross$id = NULL
  expect_identical(panel_select(y ~ x, d ~ w, cross, NULL, 'linear', control = control,
    init = 'pooled')$structure, 'cross_sectional')
  expect_error(panel_select(y ~ x, NULL, data, 'id', 'gaussian'),
    'requires a selection equation; selection = NULL fits only probit, poisson or pln')
})

test_that('subset rows equal a fit on those rows, with lm semantics for bases and padding', {
  data = interface_data(panels = 24L)
  control = interface_control()
  rows = data$id <= 20
  subset = panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', control = control,
    subset = id <= 20, init = 'pooled')
  reference = panel_select(y ~ x, d ~ x + w, data[rows, ], 'id', 'poisson', control = control,
    init = 'pooled')
  expect_same_fit(subset, reference)
  expect_identical(nrow(subset$data), nrow(data))
  # Data-dependent bases use all rows, as lm(subset =) does.
  basis = panel_select(y ~ poly(x, 2), d ~ w, data, 'id', 'poisson', control = control,
    subset = rows, init = 'pooled')
  expect_identical(attr(basis$prepared$design$outcome$terms, 'predvars'),
    attr(attr(stats::model.frame(y ~ poly(x, 2), data), 'terms'), 'predvars'))
  # Row 4 is inside the subset and row 90 (panel 23) outside it.
  gappy = data; gappy$x[c(4L, 90L)] = NA
  expect_no_error(panel_select(y ~ x, d ~ x + w, gappy, 'id', 'poisson', control = control,
    subset = id <= 20 & seq_len(nrow(gappy)) != 4L, init = 'pooled'))
  excluded = panel_select(y ~ x, d ~ x + w, gappy, 'id', 'poisson', control = control,
    subset = id <= 20, na_action = 'exclude', init = 'pooled')
  residual = residuals(excluded)
  expect_length(residual, sum(rows))
  expect_identical(names(residual), rownames(data)[rows])
  expect_true(is.na(residual[['4']]))
  expect_length(predict(excluded), sum(rows))
  missing = c(NA, rep(TRUE, nrow(data) - 1L))
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', control = control,
    subset = missing, init = 'pooled'), 'row\\(s\\): 1 \\(na_action')
  omitted = panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', control = control,
    subset = missing, na_action = 'omit', init = 'pooled')
  expect_identical(unname(unclass(omitted$prepared$na.action)), 1L)
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', subset = c(1, -2)),
    'positive or negative')
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', subset = c(1, 1)),
    'must not repeat')
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', subset = c(TRUE, FALSE)),
    'one value per data row')
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', subset = 'none'),
    'not in data')
})

test_that('specification messages name the columns, parameters and equations involved', {
  data = interface_data()
  data$x2 = 2*data$x
  control = interface_control()
  expect_error(panel_select(y ~ x + x2, d ~ x + w, data, 'id', 'poisson', control = control,
    init = 'pooled'), 'not identified by the model matrix: outcome.x2 is aliased', fixed = TRUE)
  expect_error(panel_select(y ~ x + x2, d ~ x + w, data, 'id', 'poisson', control = control,
    init = 'pooled'), 'fixed = c(`outcome.x2` = 0)', fixed = TRUE)
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', start = c(foo = 1)),
    'Unknown start parameter(s): foo; valid names: selection.(Intercept)', fixed = TRUE)
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson', start = c(cor_re = .2),
    fixed = c(cor_re = 0)), 'conflicts with a fixed parameter: cor_re')
  expect_error(panel_select(y ~ x, ~ x, data, 'id', 'poisson'), 'selection equation must be a two-sided')
  expect_error(panel_select(y ~ x, d ~ x + w, data, 'id', 'poisson',
    fixed = c(sd_re_selection = 0, cor_re = .3)), 'scale (sd_re_selection = 0)', fixed = TRUE)
  expect_error(panel_control(foo = 1), 'Unknown panel_control argument(s): foo', fixed = TRUE)
  expect_error(panel_control('agh', 10, list(), 'gh', 20, list(), 50, 1e-9, TRUE, 'optim', 'BFGS',
    500, 1e-8, 1e7, NULL, 0, .1, 'extra'), 'uniquely named')
})

test_that('control node defaults follow the route', {
  expect_identical(panel_control(nodes = NULL), panel_control())
  expect_identical(panel_control()$nodes, 10L)
  expect_identical(panel_control(quadrature = 'gl')$nodes, 20L)
  expect_identical(panel_control(observation_quadrature = 'gl')$observation_nodes, 20L)
  expect_identical(panel_control(observation_nodes = NULL)$observation_nodes, 20L)
  # Fits stored before init existed keep pooled starts when updated.
  fit = panel_select(y ~ x, d ~ x + w, interface_data(), 'id', 'poisson',
    control = interface_control(maxit = 0L), init = 'zero')
  fit$control$init = NULL
  expect_identical(update(fit, evaluate = FALSE)$init, 'pooled')
})
