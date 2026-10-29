# Observed information under a rescaled covariate. With x replaced by k*x and
# its coefficients divided by k, the likelihood is unchanged, so the exact
# information satisfies H_k = D H_1 D with D = k on the x coordinates.
information_scale_data = function() {
  local_test_seed(15)
  panels = 40L; id = rep(seq_len(panels), each = 3L); n = length(id)
  u = stats::rnorm(panels)[id]; v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)[id]
  x = stats::rnorm(n); z = stats::rnorm(n)
  d = as.integer(.2 + .5*x + .7*z + .6*u + stats::rnorm(n) > 0)
  y = as.integer(.3 + .4*x + .5*v + stats::rnorm(n) > 0)
  y[d == 0] = NA
  data.frame(id = id, x = x, z = z, d = d, y = y)
}

information_scale_case = function(k, control) {
  data = information_scale_data(); data$x = k*data$x
  prepared = prepare_data(y ~ x, d ~ x + z, data, 'id', 'probit')
  info = parameter_layout(prepared)
  fixed = resolve_fixed(c(cor_error = 0), info, prepared)
  info$fixed = info$name %in% names(fixed)
  base = c('selection.(Intercept)' = .2, selection.x = .5, selection.z = .7,
    'outcome.(Intercept)' = .3, outcome.x = .4, sd_re_selection = .6, sd_re_outcome = .5,
    cor_re = .4, cor_error = 0)
  theta = base[info$name]
  theta[c('selection.x', 'outcome.x')] = theta[c('selection.x', 'outcome.x')]/k
  evaluator = fit_evaluator(theta, info, prepared, control, fixed)
  eta = parameter_to_optimizer(theta, info)
  list(information = information_matrix(eta, evaluator), eta = eta, evaluator = evaluator)
}

expect_information_equivariance = function(control, scales) {
  unit = information_scale_case(1, control)$information
  for (k in scales) {
    scaled = information_scale_case(k, control)$information
    d = ifelse(rownames(unit) %in% c('selection.x', 'outcome.x'), k, 1)
    expected = unit*outer(d, d)
    expect_true(all(is.finite(scaled)))
    expect_lt(max(abs(c(scaled) - c(expected))/(abs(c(expected)) + 1e-8*max(abs(expected)))), 1e-8)
    expect_true(all(is.finite(positive_inverse(scaled))))
  }
}

test_that('observed information is equivariant to covariate units', {
  expect_information_equivariance(panel_control(nodes = 6L), c(1e3, 1e4))
})

test_that('checked split information is equivariant to covariate units', {
  # Split score-accuracy budgets are absolute, so much larger covariate units
  # fail the integration gate itself; that limitation is separate from the step.
  expect_information_equivariance(panel_control(quadrature = 'gl', nodes = 12L), 10)
})

test_that('observed information agrees with an independent Richardson difference', {
  case = information_scale_case(1e3, panel_control(nodes = 6L))
  gradient = function(eta) -colSums(case$evaluator$evaluate(eta, trial = FALSE)$optimizer_score)
  # Steps proportional to each coordinate's OPG standard error, independent of
  # the package's step rule.
  score = case$evaluator$evaluate(case$eta, trial = FALSE)$optimizer_score
  step = 1e-2*sqrt(diag(solve(crossprod(score))))
  central = function(h) vapply(seq_along(case$eta), function(j) {
    plus = minus = case$eta; plus[j] = plus[j] + h[j]; minus[j] = minus[j] - h[j]
    (gradient(plus) - gradient(minus))/(2*h[j])
  }, numeric(length(case$eta)))
  reference = (4*central(step/2) - central(step))/3
  reference = (reference + t(reference))/2
  expect_lt(max(abs(c(case$information) - c(reference))/max(abs(reference))), 1e-5)
})

test_that('checked information matches native near-zero stencils within roundoff', {
  # This endpoint exposed an exact-serialization lookup bug: optimHess's
  # sequential coordinate restoration can alter an untouched coordinate by ulps.
  eta = c(beta = -1.05338667132271e-7, rho = -2.48689601545857e-7)
  A = matrix(c(3, .4, .4, 2), 2L, dimnames = list(names(eta), names(eta)))
  evaluator = function(reference_factor = 1) list(
    evaluate = function(x, trial = TRUE) list(loglik = -.5*drop(crossprod(x, A %*% x)),
      optimizer_score = matrix(-drop(A %*% x), 1L, dimnames = list(NULL, names(x))),
      reference_optimizer_score = matrix(-reference_factor*drop(A %*% x), 1L,
        dimnames = list(NULL, names(x)))),
    information_gate_tolerance = 1e-8, information_steps = function(eta) rep(1e-3, 2L))
  H = information_matrix(eta, evaluator())
  expect_equal(as.vector(H), as.vector(A), tolerance = 1e-12)
  expect_identical(dim(H), dim(A))
  accuracy = attr(H, 'checked_accuracy')
  expect_identical(accuracy$status, 'PASS')
  expect_lte(accuracy$stencil_coordinate_max_scaled_roundoff, 16*.Machine$double.eps)
  expect_lte(accuracy$scaled_order_change, 1e-12)
  expect_lte(accuracy$scaled_max_change, 1e-12)
  for (match in accuracy$stencil_coordinate_matches) {
    expect_equal(match$actual - match$requested, match$difference, tolerance = 0)
    expect_lte(max(abs(match$difference)/pmax(1, abs(match$actual), abs(match$requested))),
      16*.Machine$double.eps)
  }
  # Step stability alone cannot certify information: a persistent lower-order
  # derivative discrepancy must leave covariance unavailable.
  expect_warning({unresolved = information_matrix(eta, evaluator(1.01)); unresolved}, 'unavailable')
  expect_true(all(is.na(unresolved)))
})

test_that('checked observed covariance rejects unstable inverses and preserves non-PD diagnostics', {
  # Entry changes are standardized by the diagonals, so a curvature
  # change in one weak coordinate fails the entry check itself.
  expect_error(information_accuracy_gate(diag(c(1, 1e-8)), diag(c(1, 1e-12)),
    diag(c(1, 1e-8)), diag(c(1, 1e-12)), .002),
    'step/order refinement unresolved', class = 'panelselect_integration_error')
  correlated = function(rho) matrix(c(1, rho, rho, 1), 2L)
  native = correlated(1 - 1e-8)
  finer = correlated(1 - 1e-12)
  # Standardized entry changes below 1e-8 hide a 10,000-fold variance change
  # near singularity. Positive definiteness and derivative agreement alone are
  # insufficient.
  expect_lt(max(abs(native - finer)), 1e-7)
  expect_error(information_accuracy_gate(native, finer, native, finer, .002),
    'covariance refinement unresolved', class = 'panelselect_integration_error')
  indefinite = correlated(1 + 1e-12)
  expect_error(information_accuracy_gate(native, native, indefinite, native, .002),
    'not positive definite', class = 'panelselect_integration_error')
  # Retain stable information at arbitrary nonmaximum vectors; it is useful
  # derivative evidence but must not be described as an available covariance.
  non_pd = diag(c(1, -1))
  value = information_accuracy_gate(non_pd, non_pd, non_pd, non_pd, .002)
  expect_identical(value$status, 'DERIVATIVE_REFINEMENT_PASS_COVARIANCE_UNAVAILABLE')
  expect_identical(value$observed_covariance_refinement$status, 'UNAVAILABLE_NATIVE_NOT_PD')
  expect_identical(value$canonical_information, non_pd)
  regular = diag(c(2, 3)); close = diag(c(2.0001, 3.0001))
  value = information_accuracy_gate(regular, close, regular, close, .002)
  expect_identical(value$status, 'PASS')
  expect_identical(value$observed_covariance_refinement$status, 'PASS')
  expect_length(value$observed_covariance_refinement$comparisons, 3L)
})
