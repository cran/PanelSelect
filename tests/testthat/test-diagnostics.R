test_that("quadrature checks refine PLN rules independently without mutating the fit", {
  data = data.frame(id = rep(1:2, each = 3), y = c(0, 1, 2, 1, 0, 3))
  fit = panel_select(y ~ 1, NULL, data, "id", "pln",
    fixed = c("outcome.(Intercept)" = .2, sd_re_outcome = .3, sd_error = .2),
    control = panel_control(nodes = 3L, observation_nodes = 5L))
  saved = fit
  check = check_quadrature(fit, nodes = 5L, observation_nodes = 9L)
  expect_equal(check$comparison$rule, c("person_effect", "observation", "both"))
  expect_false(any(check$comparison$integration_failed))
  expect_null(check$refit)
  expect_identical(fit, saved)
  expect_error(check_quadrature(fit, nodes = 2L), "not be smaller")
  expect_error(check_quadrature(fit, nodes = 3L, observation_nodes = 5L), "must increase")
  refined = check_quadrature(fit, nodes = 5L, observation_nodes = 9L, refit = TRUE)$refit
  expect_equal(refined$control$nodes, 5L)
  expect_equal(refined$control$observation_nodes, 9L)
  expect_equal(coef(refined), coef(fit))
  # All parameters are fixed, so the refit log likelihood moves by exactly the
  # joint-rule change that the comparison reports.
  expect_equal(as.numeric(logLik(refined)) - as.numeric(logLik(fit)),
    check$comparison$loglik_change[check$comparison$rule == "both"], tolerance = 1e-12)
})

test_that("quadrature checks require an effective integration order increase", {
  data = data.frame(id = rep(1:3, each = 3), y = rep(c(0, 1, 2), 3))
  fixed = c("outcome.(Intercept)" = .2, sd_re_outcome = 0)
  control = panel_control(quadrature = "gh", nodes = c(3L, 5L), observation_nodes = 5L)
  direct = panel_select(y ~ 1, NULL, data, "id", "poisson", fixed = fixed, control = control)
  expect_error(check_quadrature(direct), "no numerical integral")
  fixed["sd_re_outcome"] = .3
  one = panel_select(y ~ 1, NULL, data, "id", "poisson", fixed = fixed, control = control)
  # A standalone model's only effect is the outcome effect: element 2 of nodes.
  expect_identical(check_quadrature(one, nodes = c(3L, 7L))$comparison$rule, "person_effect")
  expect_error(check_quadrature(one, nodes = c(3, 5)), "must increase")
  expect_error(check_quadrature(one, nodes = c(5L, 5L)), "effective quadrature order")
  expect_error(check_quadrature(one, nodes = c(7L, 5L)), "effective quadrature order")

  fixed["sd_re_outcome"] = 0
  inner = panel_select(y ~ 1, NULL, data, "id", "pln", fixed = c(fixed, sd_error = .3),
    control = control)
  expect_identical(check_quadrature(inner)$comparison$rule, "observation")
  fixed["sd_re_outcome"] = .3
  both = panel_select(y ~ 1, NULL, data, "id", "pln", fixed = c(fixed, sd_error = .3),
    control = control)
  expect_identical(check_quadrature(both, nodes = c(3L, 7L))$comparison$rule,
    c("person_effect", "observation", "both"))
})

test_that("quadrature checks follow the orders of independent factors", {
  data = data.frame(id = rep(1:3, each = 3), y = rep(c(0, 1, 2), 3), d = 1)
  fixed = c("selection.(Intercept)" = .1, "outcome.(Intercept)" = .2,
    sd_re_selection = 0, sd_re_outcome = .3, cor_re = 0)
  fit = panel_select(y ~ 1, d ~ 1, data, "id", "poisson", fixed = fixed,
    control = panel_control(quadrature = "gh", nodes = c(3L, 5L)))
  # Only the outcome factor integrates, using the second supplied order.
  expect_error(check_quadrature(fit, nodes = c(7L, 5L)), "effective quadrature order")
  expect_identical(check_quadrature(fit, nodes = c(3L, 7L))$comparison$rule, "person_effect")
})

test_that("unobserved PLN outcomes do not create an inner refinement", {
  data = data.frame(id = rep(1:3, each = 3), y = NA_real_, d = 0)
  fixed = c("selection.(Intercept)" = .1, "outcome.(Intercept)" = .2,
    sd_re_selection = .3, sd_re_outcome = .4, sd_error = .5, cor_re = .2, cor_error = .1)
  fit = panel_select(y ~ 1, d ~ 1, data, "id", "pln", fixed = fixed,
    control = panel_control(quadrature = "gh", nodes = c(3L, 5L), observation_nodes = 5L))
  expect_error(check_quadrature(fit, nodes = c(3L, 5L), observation_nodes = 7L), "must increase")
  expect_identical(check_quadrature(fit, nodes = c(5L, 5L))$comparison$rule, "person_effect")
})
