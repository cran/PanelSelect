# Historical fits use an earlier vocabulary. A current fit rewritten in that
# vocabulary must behave exactly as the current fit in every method.

# Person effects correlated 0.4 and errors correlated 0.3: interior values.
vocabulary_data = function(seed = 1401, panels = 24L, waves = 4L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  e = stats::rnorm(n); epsilon = .3*e + sqrt(1 - .3^2)*stats::rnorm(n)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + e > 0)
  data$y = ifelse(data$d == 1, .5 + .4*data$x + .6*v[data$id] + .5*epsilon, NA)
  data$count = stats::rpois(n, exp(.1 + .3*data$x + .5*u[data$id]))
  data
}

# Staged integration diagnostics record elapsed seconds, which vary by run.
without_timings = function(x) {
  if (!is.list(x)) return(x)
  x[names(x) %in% 'elapsed_seconds'] = NULL
  for (i in seq_along(x)) if (is.list(x[[i]])) x[[i]] = without_timings(x[[i]])
  x
}

expect_same_methods = function(old, fit, data, history = TRUE) {
  rows = data[1:4, ]
  normalized = current_fit(old); normalized$call = fit$call
  expect_identical(normalized, fit)
  expect_identical(coef(old), coef(fit))
  expect_identical(vcov(old), vcov(fit))
  expect_identical(vcov(old, 'observed'), vcov(fit, 'observed'))
  expect_identical(confint(old), confint(fit))
  expect_identical(utils::capture.output(print(old)), utils::capture.output(print(fit)))
  expect_identical(utils::capture.output(print(summary(old))),
    utils::capture.output(print(summary(fit))))
  expect_identical(predict(old, rows), predict(fit, rows))
  expect_identical(predict(old, rows, se.fit = TRUE), predict(fit, rows, se.fit = TRUE))
  expect_identical(residuals(old), residuals(fit))
  informed = add_information(old); informed$call = fit$call
  expect_identical(informed, add_information(fit))
  expect_identical(check_stationarity(old), check_stationarity(fit))
  expect_identical(without_timings(check_quadrature(old)), without_timings(check_quadrature(fit)))
  expect_identical(update(old, evaluate = FALSE), update(fit, evaluate = FALSE))
  if (history) {
    expect_identical(predict(old, rows, conditioning = 'history', history = data),
      predict(fit, rows, conditioning = 'history', history = data))
    expect_identical(random_effects(old), random_effects(fit))
  }
}

test_that('a Gauss-Hermite fit saved in the earlier vocabulary behaves as the current fit', {
  data = vocabulary_data()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'linear',
    control = panel_control(quadrature = 'agh', nodes = 4L),
    init = 'pooled')
  old = old_vocabulary(fit)
  expect_true(all(c('inner_quadrature', 'split', 'reduction', 'stationarity_tol', 'vcov') %in%
    names(old$control)))
  expect_identical(old$control$quadrature, 'adaptive')
  expect_true('sd_selection' %in% names(old$coefficients))
  expect_true('loglik_by_panel' %in% names(old))
  expect_same_methods(old, fit, data)
  # A historical check named its comparison method and observation columns.
  check = check_quadrature(fit)
  saved = check; saved$method = 'refine'; saved$compare = NULL
  names(saved$comparison) = sub('^observation_', 'inner_', names(saved$comparison))
  saved$comparison$inner_quadrature[!is.na(saved$comparison$inner_quadrature)] = 'fixed'
  expect_identical(utils::capture.output(print(saved)), utils::capture.output(print(check)))
  # The earlier record names its arguments in the earlier order.
  record = update(old_vocabulary(fit, formula_record = TRUE), evaluate = TRUE)$call
  expect_identical(names(record)[2:6], c('outcome', 'selection', 'data', 'id', 'family'))
  # An earlier record that relied on the gaussian default and gave the fit
  # settings inside panel_control() is written as a current call.
  earlier = old_vocabulary(fit)
  earlier$call$family = NULL
  earlier$call$control = quote(panel_control(nodes = 4L, init = 'pooled', vcov = 'opg'))
  record = update(earlier)$call
  expect_identical(record$family, 'gaussian')
  expect_identical(record$control, quote(panel_control(nodes = 4L)))
  expect_identical(c(record$init, record$vcov), c('pooled', 'opg'))
})

test_that('a staged GL fit and its retained mesh saved in the earlier vocabulary behave as the current fit', {
  data = vocabulary_data()
  fit = panel_select(count ~ x, NULL, data, 'id', 'poisson',
    control = staged_test_budget(panel_control(quadrature = 'gl', nodes = 12L)))
  expect_false(is.null(fit$control$integration_mesh))
  old = old_vocabulary(fit)
  expect_identical(old$control$quadrature, 'split')
  expect_true(all(c('loglik_tol', 'max_intervals') %in% names(old$control$split)))
  expect_identical(old$control$integration_mesh$binding$outer, 'split')
  expect_true('reduction' %in% names(old$control$integration_mesh$binding))
  expect_identical(old$control$integration_mesh$outer[[1L]]$method, 'split')
  expect_same_methods(old, fit, data, history = FALSE)
})

test_that('a fixed-zero standalone fit saved in the earlier vocabulary keeps its restriction', {
  data = vocabulary_data()
  fit = panel_select(count ~ x, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0))
  old = old_vocabulary(fit)
  expect_identical(old$fixed, c(sd_outcome = 0))
  expect_identical(coef(update(old)), coef(update(fit)))
  expect_identical(update(old, evaluate = FALSE)$fixed, c(sd_re_outcome = 0))
})
