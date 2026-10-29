test_that("selection-only panel splits preserve computed component diagnostics", {
  fixture = core_fixture("selection_pln")
  prepared = fixture$prepared
  prepared$group = c(1L, 2L, 2L); prepared$group_ptr = c(0L, 1L, 3L)
  prepared$id_values = c("selection-only", "observed"); prepared$n_panels = 2L
  prepared$d = c(0, 1, 1)
  calls = list()
  local_mocked_bindings(integrate_prepared = function(theta, prepared, control, fixed, score = TRUE) {
    kind = if (prepared$has_selection) "joint" else prepared$family
    diagnostics = list(method = "agh", component = kind,
      mode_residual = rep(1e-12, prepared$n_panels), nested = list(nodes = control$nodes))
    calls[[length(calls)+1L]] <<- diagnostics
    value = if (kind == "joint") -5 else -3
    list(loglik = rep(value, prepared$n_panels),
         score = matrix(1, prepared$n_panels, length(theta), dimnames = list(NULL, names(theta))),
         diagnostics = diagnostics)
  }, .package = "PanelSelect")
  answer = evaluate_model(fixture$theta, prepared, core_control("agh"))
  expect_length(calls, 2L)
  expect_identical(answer$loglik, c(-3, -5))
  expected = matrix(1, 2L, length(fixture$theta))
  expected[1L, ] = as.numeric(startsWith(names(fixture$theta), 'selection.') |
                               names(fixture$theta) == 'sd_re_selection')
  expect_equal(unname(answer$score), expected)
  expect_identical(answer$diagnostics$selection_only, calls[[1L]])
  expect_identical(answer$diagnostics$observed, calls[[2L]])
  expect_identical(answer$diagnostics$selection_only_panels, "selection-only")
  expect_identical(answer$diagnostics$observed_panels, "observed")

  calls = list(); fixed = c(cor_re = 0, cor_error = 0)
  fixture$theta[names(fixed)] = fixed
  answer = evaluate_model(fixture$theta, prepared, core_control("agh"), fixed)
  expect_length(calls, 3L)
  expect_true(answer$diagnostics$observed$factorized)
  expect_identical(answer$diagnostics$observed$selection, calls[[2L]])
  expect_identical(answer$diagnostics$observed$outcome, calls[[3L]])

  calls = list(); prepared$d[] = 0
  answer = evaluate_model(fixture$theta, prepared, core_control("agh"), fixed)
  expect_length(calls, 1L)
  expect_identical(answer$diagnostics$selection_only, calls[[1L]])
  expect_null(answer$diagnostics$observed)
  expect_identical(answer$diagnostics$observed_panels, character())
})

test_that("observed information error text is retained without replacing warnings or covariance", {
  perturbations = 0L
  local_mocked_bindings(fit_evaluator = function(template, info, prepared, control, fixed) {
    list(evaluate = function(eta, trial = TRUE) {
      if (any(abs(eta) > 1e-8)) {
        perturbations <<- perturbations+1L
        integration_error("fixture perturbed mode failure")
      }
      list(loglik = rep(-1, prepared$n_panels),
        score = matrix(0, prepared$n_panels, nrow(info)),
        optimizer_score = matrix(0, prepared$n_panels, sum(!info$fixed)),
        theta = template, valid = TRUE, diagnostics = list(method = "gh"))
    }, rejected = function() 0L, information_steps = function(eta) rep(1e-3, sum(!info$fixed)))
  }, .package = "PanelSelect")
  warnings = character()
  fit = withCallingHandlers(panel_select(y ~ 1, NULL, data.frame(id = 1:3, y = 0:2), "id", "poisson",
    start = c("outcome.(Intercept)" = 0), fixed = c(sd_re_outcome = 0),
    control = panel_control(maxit = 0L), vcov = "observed"), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_identical(fit$optimizer$information_error, "fixture perturbed mode failure")
  expect_identical(warnings, c("Observed information is unavailable: fixture perturbed mode failure",
    "The requested covariance estimate is unavailable; inspect the fit diagnostics"))
  expect_identical(perturbations, 1L)
  expect_true(fit$optimizer$information_computed)
  expect_true(all(is.na(fit$optimizer$information)))
  expect_true(is.na(fit$covariance$observed["outcome.(Intercept)", "outcome.(Intercept)"]))
})

test_that("successful information diagnostics add no evaluations or matrix changes", {
  calls = 0L; message = NULL
  evaluator = list(evaluate = function(eta, trial = TRUE) {
    calls <<- calls+1L
    list(loglik = -sum(eta^2)/2, optimizer_score = matrix(-eta, 1L))
  }, information_steps = function(eta) rep(1e-3, 2L))
  eta = c(a = .2, b = -.1)
  baseline = information_matrix(eta, evaluator)
  baseline_calls = calls; calls = 0L
  observed = information_matrix(eta, evaluator, on_error = function(error) message <<- error)
  expect_identical(observed, baseline)
  expect_identical(calls, baseline_calls)
  expect_null(message)
})
