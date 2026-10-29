# Common mistakes stop with a message that names the argument involved and,
# where one exists, the fix. Base R messages (unused arguments, missing
# variables, ambiguous partial names) are accepted where they already do so.

mistakes_data = function(seed = 1701, panels = 20L, waves = 4L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n), w = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + stats::rnorm(n) > 0)
  data$y = ifelse(data$d == 1, .5 + .4*data$x + .6*v[data$id] + stats::rnorm(n), NA)
  data$count = stats::rpois(n, exp(.2 + .3*data$x + .5*u[data$id]))
  data$b = as.integer(data$x + stats::rnorm(n) > 0)
  data
}

test_that('fitting-call mistakes name the argument and the fix', {
  data = mistakes_data()
  expect_error(panel_select(y ~ x, data, 'id', 'linear', d ~ w),
    'takes \\(outcome, selection, data, id, family\\)')
  expect_error(panel_select(y ~ x, d ~ w, data, 'id'), 'family is required')
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'binomial'),
    'family = "binomial" is not available; use "linear"', fixed = TRUE)
  expect_error(panel_select(b ~ x, NULL, data, 'id', 'logit'), 'family = "logit" is not available',
    fixed = TRUE)
  expect_error(panel_select(b ~ x, NULL, data, 'id', binomial()), 'a family object such as binomial')
  expect_error(panel_select(y ~ x, d ~ w, data, 'person', 'linear'), 'id must name one column in data')
  expect_error(panel_select(y ~ x, d ~ w, data, data$id, 'linear'), 'id must name one column in data')
  expect_error(panel_select(y ~ x, ~ w, data, 'id', 'linear'), 'selection equation must be a two-sided')
  expect_error(panel_select(z ~ x, d ~ w, data, 'id', 'linear'), "object 'z' not found")
  expect_error(panel_select(d ~ x, NULL, data, 'id', 'linear'), 'requires a selection equation')
  expect_error(panel_select(count ~ x, NULL, data, 'id', 'probit'), 'Binary outcomes must be zero or one')
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'poisson'), 'Count outcomes must be nonnegative integers')
  expect_error(panel_select(count ~ x, NULL, data, 'id', 'poisson', maxit = 100),
    'unused argument (maxit = 100)', fixed = TRUE)
  expect_error(panel_select(count ~ x, NULL, data, 'id', 'poisson', init = 'probit'), 'init must be one of')
  expect_error(panel_select(count ~ x, NULL, transform(data, x = replace(x, 1L, NA)), 'id', 'poisson'),
    'Missing or nonfinite required data in row(s): 1', fixed = TRUE)
})

test_that('parameter mistakes name the parameter and the valid names', {
  data = mistakes_data()
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'linear', start = c(foo = 1)),
    'Unknown start parameter(s): foo; valid names: selection.(Intercept)', fixed = TRUE)
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'linear', start = c(sd_selection = 1)),
    'sd_selection is now sd_re_selection (in start)', fixed = TRUE)
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'linear', fixed = c(cor_re = 1)),
    'fixed correlations must lie strictly between -1 and 1')
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'linear', fixed = c(sd_re_outcome = -1)),
    'fixed standard deviations must be nonnegative')
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', 'linear', start = c(cor_re = 1.5)),
    'start correlations must lie strictly between -1 and 1')
  # Cross-sectional models have no person-effect SDs.
  expect_error(panel_select(y ~ x, d ~ w, data, NULL, 'linear', start = c(sd_re_outcome = 1)),
    'Unknown start parameter(s): sd_re_outcome', fixed = TRUE)
})

test_that('control mistakes name the option and the rule', {
  expect_error(panel_control(nodes = 0), 'nodes must contain valid integer values')
  expect_error(panel_control(quadrature = 'gl', quadrature_options = list(bounds = 8)),
    'Unknown quadrature_options entry for the gl rule: bounds', fixed = TRUE)
  expect_error(panel_control(quadrature_options = list(bound = 8)),
    "quadrature_options require quadrature = 'gl'", fixed = TRUE)
  expect_error(panel_control(quadrature = 'adaptive'), 'quadrature = "adaptive" is now "agh"', fixed = TRUE)
  expect_error(panel_control(vcov = 'observed'), 'vcov is now an argument of panel_select()', fixed = TRUE)
  # Prefixes shared by several arguments are ambiguous.
  expect_error(panel_control(observation = 'gl'), 'matches multiple formal arguments')
  expect_error(panel_control(quad = 'gl'), 'matches multiple formal arguments')
})

test_that('method mistakes name the argument', {
  data = mistakes_data()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'linear', control = panel_control(nodes = 4L),
    init = 'pooled')
  rows = data[1:2, ]
  expect_error(predict(fit, rows, type = 'probability'), 'should be one of')
  expect_error(predict(fit, rows, type = 'selection'), 'type = "selection" is now "selection_probability"',
    fixed = TRUE)
  expect_error(summary(fit, vcov = 'robust'), 'should be one of')
  expect_error(confint(fit, 'sd_selection'), 'sd_selection is now sd_re_selection (in parm)', fixed = TRUE)
  expect_error(confint(fit, 'foo'), 'Unknown parameter in parm: foo', fixed = TRUE)
  expect_error(update(fit, maxit = 10), 'Update overrides must name panel_select arguments')
  expect_error(check_quadrature(fit, nodes = 2), 'Refined orders must not be smaller than fitted orders')
  expect_error(check_quadrature(fit, method = 'refine'), 'method is now compare')
})
