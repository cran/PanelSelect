serialization_fixture = function() {
  hidden_data = data.frame(id = rep(seq_len(6), each = 2),
    y = c(1, 2, 0, 3, 1, 0, 2, 1, 4, 0, 1, 3),
    x = c(-.8, .2, -.4, .7, .1, -.2, .5, -.3, .9, -.6, .3, .4))
  hidden_formula = y ~ x
  panel_select(hidden_formula, NULL, hidden_data, "id", "poisson", fixed = c(sd_re_outcome = 0))
}

test_that("serialized models retain post-estimation and self-contained updates", {
  fit = serialization_fixture()
  path = tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(fit, path)
  fresh = new.env(parent = baseenv())
  fresh$fit = readRDS(path)
  expect_equal(eval(quote(stats::coef(fit)), fresh), coef(fit))
  expect_equal(eval(quote(stats::vcov(fit)), fresh), vcov(fit))
  expect_equal(eval(quote(summary(fit)$coefficients), fresh), summary(fit)$coefficients)
  expect_equal(eval(quote(stats::predict(fit)), fresh), predict(fit))
  expect_equal(eval(quote(stats::residuals(fit)), fresh), residuals(fit))
  expect_equal(eval(quote(stats::model.matrix(fit)), fresh), model.matrix(fit))
  updated = eval(quote(stats::update(fit)), fresh)
  expect_equal(coef(updated), coef(fit), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(updated)), as.numeric(logLik(fit)), tolerance = 1e-8)
  call = eval(quote(stats::update(fit, evaluate = FALSE)), fresh)
  expect_identical(call[[1L]], quote(PanelSelect::panel_select))
  expect_identical(call$data, fit$data)
  expect_equal(coef(eval(call, fresh)), coef(fit), tolerance = 1e-6)
})

test_that("updates resolve overrides and remove only obsolete default parameters", {
  fit = serialization_fixture()
  fresh = new.env(parent = baseenv())
  fresh$fit = fit
  fresh$replacement = transform(fit$data, x = x + .1)
  call = eval(quote(stats::update(fit, . ~ 1, data = replacement,
    fixed = c(sd_re_outcome = .2), evaluate = FALSE)), fresh)
  expect_identical(call$data, fresh$replacement)
  expect_false("outcome.x" %in% names(call$start))
  expect_equal(call$start[["sd_re_outcome"]], .2)
  expect_s3_class(eval(call, fresh), "panelselect")
  released = update(fit, fixed = NULL, evaluate = FALSE)
  expect_false("sd_re_outcome" %in% names(released$start))
  expect_error(update(fit, fixed = NULL, start = c(sd_re_outcome = 0)), "must be positive")
  expect_error(update(fit, . ~ 1, start = c(outcome.x = 0)), "Unknown start")
  expect_error(update(fit, . ~ 1, fixed = c(outcome.x = 0)), "Unknown fixed")
  changed = update(fit, family = "pln", evaluate = FALSE)
  expect_false("sd_error" %in% names(changed$start))
  expect_identical(changed$family, "pln")
  pln = panel_select(y ~ x, NULL, fit$data, "id", "pln", fixed = c("outcome.(Intercept)" = .3, outcome.x = .2, sd_re_outcome = .4, sd_error = .2))
  reduced = update(pln, family = "poisson", evaluate = FALSE)
  expect_false("sd_error" %in% names(reduced$start))
  expect_false("sd_error" %in% names(reduced$fixed))
  expect_s3_class(eval(reduced, fresh), "panelselect")
})

test_that("updates after serialization do not need vanished call locals", {
  # A call recorded inside another function can name locals that are gone
  # once the fit is saved.
  fit = serialization_fixture()
  fit$call = quote(panel_select(local_formula, NULL, vanished_data, id_column,
                                 family = local_family))
  path = tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(fit, path)
  fresh = new.env(parent = baseenv())
  fresh$fit = readRDS(path)
  updated = eval(quote(stats::update(fit)), fresh)
  expect_s3_class(updated, "panelselect")
  expect_equal(coef(updated), coef(fit), tolerance = 1e-6)
  expect_equal(predict(updated), predict(fit), tolerance = 1e-6)
})
