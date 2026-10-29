# Sandwich covariance: an independent GLM cluster-robust reference and an
# exact assembly identity from the stored information and panel scores.

test_that('panel sandwich covariance equals the clustered GLM sandwich at sd_re_outcome = 0', {
  local_test_seed(8101)
  data = data.frame(id = rep(1:12, each = 3), x = stats::rnorm(36))
  data$y = stats::rpois(36, exp(.3 + .5*data$x + stats::rnorm(12)[data$id]*.4))
  fit = panel_select(y ~ x, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0), vcov = 'sandwich')
  beta = coef(fit)[c('outcome.(Intercept)', 'outcome.x')]
  X = cbind(1, data$x)
  mu = drop(exp(X %*% beta))
  information = crossprod(X, X*mu)
  scores = rowsum(X*(data$y - mu), data$id)
  inverse = solve(information)
  expected = inverse %*% crossprod(scores) %*% inverse
  actual = vcov(fit, type = 'sandwich')[names(beta), names(beta)]
  expect_equal(unname(actual), unname(expected), tolerance = 1e-6)
  expect_identical(fit$optimizer$status, 'converged')
})

test_that('sandwich covariance assembles H^-1 (sum of panel score products) H^-1 on optimizer scale', {
  local_test_seed(8102)
  data = data.frame(id = rep(1:15, each = 3), x = stats::rnorm(45), w = stats::rnorm(45))
  u = stats::rnorm(15)[data$id]
  data$d = as.integer(.2 + .6*data$w + .7*u + stats::rnorm(45) > 0)
  data$y = ifelse(data$d == 1, .4 + .5*data$x + .5*u + stats::rnorm(45), NA)
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'gaussian', fixed = c(cor_error = 0),
    control = panel_control(nodes = 6L), vcov = 'sandwich')
  theta = coef(fit); free = !fit$parameter_info$fixed
  type = fit$parameter_info$type
  # Natural-to-optimizer derivatives written out: log scales and atanh correlations.
  jacobian = ifelse(type == 'sd', theta, ifelse(type == 'correlation', 1 - theta^2, 1))[free]
  scores = sweep(fit$score[, free, drop = FALSE], 2L, jacobian, `*`)
  inverse = solve(fit$optimizer$information)
  expected = (inverse %*% crossprod(scores) %*% inverse) * outer(jacobian, jacobian)
  actual = vcov(fit, type = 'sandwich')[free, free]
  expect_equal(unname(actual), unname(expected), tolerance = 1e-10)
})
