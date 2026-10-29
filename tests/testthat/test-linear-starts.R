# Linear selection data simulated from the model; the global RNG state is
# restored afterwards. `scale` multiplies the outcome.
linear_scale_data = function(scale = 1) {
  local_test_seed(14)
  panels = 60L; id = rep(seq_len(panels), each = 4L); n = length(id)
  u = stats::rnorm(panels)[id]; v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)[id]
  x = stats::rnorm(n); z = stats::rnorm(n)
  d = as.integer(.2 + .5*x + .7*z + .6*u + stats::rnorm(n) > 0)
  y = .3 + .4*x + .5*v + stats::rnorm(n)
  y[d == 0] = NA
  data.frame(person = id, x = x, z = z, d = d, y = scale*y)
}

linear_scale_fit = function(data, ...) panel_select(y ~ x,
  d ~ x + z, data, 'person', 'gaussian', fixed = c(cor_error = 0), control = panel_control(nodes = 6L, ...))

test_that('Linear scale starts follow the pooled residual SD of the selected outcomes', {
  data = linear_scale_data()
  prepared = prepare_data(y ~ x, d ~ x + z, data, 'person', 'gaussian')
  info = parameter_layout(prepared)
  fixed = resolve_fixed(c(cor_error = 0), info, prepared)
  info$fixed = info$name %in% names(fixed)
  ols = stats::lm(y ~ x, data = data[data$d == 1, ])
  s = sqrt(sum(stats::residuals(ols)^2)/ols$df.residual)
  theta = initial_parameters(prepared, info, NULL, fixed, fit_control(init = 'pooled'))
  expect_equal(unname(theta[c('sd_re_outcome', 'sd_error')]), rep(s/sqrt(2), 2), tolerance = 1e-12)
  # A supplied scale prevails and the remaining free scale carries the variance.
  supplied = initial_parameters(prepared, info, c(sd_re_outcome = .7), fixed, fit_control(init = 'pooled'))
  expect_identical(unname(supplied['sd_re_outcome']), .7)
  expect_equal(unname(supplied['sd_error']), s, tolerance = 1e-12)
})

test_that('Linear selection fits are equivariant to the outcome scale', {
  unit = linear_scale_fit(linear_scale_data(), reltol = 1e-14)
  scaled = linear_scale_fit(linear_scale_data(1000), reltol = 1e-14)
  selected = sum(linear_scale_data()$d == 1)
  # Exact identity: the density of 1000*y is the density of y divided by 1000.
  expect_lt(abs(scaled$loglik - (unit$loglik - selected*log(1000))), 1e-6)
  multiplier = ifelse(startsWith(names(coef(unit)), 'outcome.') |
    names(coef(unit)) %in% c('sd_re_outcome', 'sd_error'), 1000, 1)
  expect_equal(coef(scaled), coef(unit)*multiplier, tolerance = 1e-5)
  default = linear_scale_fit(linear_scale_data(1000))
  expect_identical(default$optimizer$status, 'converged')
  expect_lt(abs(default$loglik - (unit$loglik - selected*log(1000))), .5)
})
