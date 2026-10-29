sharing_fixture = function(family, na_action = 'fail') {
  data = data.frame(id = c(2, 1, 2, 1, 3, 3), x = c(-1, 0, .2, .5, 1, -.4),
    s = c(.1, -.2, .4, .3, -.5, .8), off = c(.1, -.1, .2, 0, .1, .15),
    d = c(1, 0, 1, 1, 0, 1), y = c(1, NA, 0, 1, NA, 0))
  if (family == 'gaussian') data$y = data$y+.3*data$x
  if (na_action != 'fail') data$x[3L] = NA_real_
  prepared = prepare_data(y ~ x+offset(off), d ~ s, data, 'id', family, na_action)
  info = parameter_layout(prepared)
  theta = stats::setNames(rep(.2, nrow(info)), info$name)
  theta[info$type == 'sd'] = .6
  theta[info$type == 'correlation'] = -.25
  V = diag(.01, length(theta)); dimnames(V) = list(names(theta), names(theta))
  structure(list(coefficients = theta, prepared = prepared, parameter_info = info,
    covariance = list(observed = V, opg = 2*V, sandwich = 3*V),
    control = fit_control(), fixed = numeric(), data = data,
    optimizer = list(information_computed = TRUE)), class = 'panelselect')
}

test_that('marginal predictions avoid joint work and resolve covariance once', {
  fit = sharing_fixture('probit')
  component_calls = design_calls = covariance_calls = 0L
  components = model_predictions; design = prediction_design; covariance = prediction_covariance
  testthat::local_mocked_bindings(
    model_predictions = function(...) { component_calls <<- component_calls+1L; components(...) },
    prediction_design = function(...) { design_calls <<- design_calls+1L; design(...) },
    prediction_covariance = function(...) { covariance_calls <<- covariance_calls+1L; covariance(...) },
    prediction_bvn = function(...) stop('An unused joint target was evaluated'))
  # Each marginal target builds its own equation's design only.
  for (type in c('selection_probability', 'response')) {
    component_calls = design_calls = covariance_calls = 0L
    result = predict(fit, fit$data, type = type, se.fit = TRUE)
    expect_identical(c(component_calls, design_calls, covariance_calls), c(1L, 1L, 1L), label = type)
    expect_true(all(is.finite(result$se.fit)), label = type)
  }
  fit$covariance[[fit$control$vcov]][] = NA_real_
  for (type in c('selection_probability', 'response'))
    expect_true(all(is.na(predict(fit, fit$data, type = type, se.fit = TRUE)$se.fit)), label = type)
  expect_true(is.na(predict(fit, fit$data, aggregate = TRUE, se.fit = TRUE)$se.fit))
})

test_that('vectorized prediction variances match scalar full-covariance forms', {
  for (n in c(1L, 4L, 9L)) for (p in c(1L, 3L, 5L)) {
    gradient = matrix(sin(seq_len(n*p))*rep(seq_len(p), each = n), n, p)
    root = matrix(cos(seq_len(p*p)), p, p)
    covariance = crossprod(root)+diag(.2, p)
    expected = vapply(seq_len(n), function(i)
      sqrt(drop(gradient[i, , drop = FALSE] %*% covariance %*% gradient[i, ])), numeric(1))
    expect_equal(prediction_se(gradient, covariance), expected, tolerance = 1e-13)
  }
  gradient = rbind(c(1, 0), c(0, 1), c(0, 0), c(NA, 1))
  covariance = matrix(c(4, NA, NA, 9), 2L)
  expect_equal(prediction_se(gradient, covariance), c(2, 3, 0, NA))
  covariance = diag(c(-1, 4))
  expect_equal(prediction_se(gradient, covariance), c(NA, 2, 0, NA))
  covariance = diag(c(-1e-14, 4))
  expect_equal(prediction_se(gradient[1:3, ], covariance), c(0, 2, 0))
})
