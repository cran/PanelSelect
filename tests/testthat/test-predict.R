prediction_fixture = function(family = "probit", selection = TRUE, adaptive = FALSE) {
  data = data.frame(id = rep(1:3, each = 2), x = c(-1, 0, .2, .5, 1, -.4),
                     off = c(.1, -.1, .2, 0, .1, .15), d = c(1, 0, 1, 1, 0, 1),
                     y = c(1, 0, 0, 1, 0, 1))
  if (family == "gaussian") data$y = data$y + .3*data$x
  if (family %in% c("poisson", "pln")) data$y = c(2, 1, 0, 3, 2, 1)
  if (selection) data$y[data$d == 0] = NA_real_
  prepared = prepare_data(y ~ x + offset(off), if (selection) d ~ x else NULL,
                            data, "id", family)
  info = parameter_layout(prepared)
  theta = stats::setNames(rep(.2, nrow(info)), info$name)
  theta[info$type == "sd"] = .6
  theta[info$type == "correlation"] = -.25
  if ("cor_error" %in% names(theta)) theta["cor_error"] = .2
  V = diag(rep(.01, length(theta))); dimnames(V) = list(names(theta), names(theta))
  structure(list(coefficients = theta, prepared = prepared, parameter_info = info,
    covariance = list(observed = V, opg = 2*V, sandwich = 3*V),
    control = fit_control(panel_control(quadrature = if (adaptive) "agh" else "gh",
                                        nodes = 6, observation_nodes = 12)),
    fixed = numeric(), data = data, optimizer = list(information_computed = TRUE)),
    class = "panelselect")
}

prediction_difference = function(fit, fun, step = 1e-5) {
  theta = fit$coefficients
  vapply(seq_along(theta), function(j) {
    high = low = fit
    high$coefficients[j] = theta[j] + step
    low$coefficients[j] = theta[j] - step
    (fun(high)-fun(low))/(2*step)
  }, numeric(length(fun(fit))))
}

test_that("population prediction gradients match independent perturbations for every model", {
  for (family in c("gaussian", "probit", "poisson", "pln"))
    for (selection in if (family == "gaussian") TRUE else c(FALSE, TRUE)) {
      fit = prediction_fixture(family, selection)
      targets = c("response", "link", if (selection) c("selection_probability", "joint", "selected_mean"))
      for (type in targets) {
        result = predict(fit, type = type, se.fit = TRUE)
        expect_identical(colnames(result$gradient), names(fit$coefficients))
        numeric = prediction_difference(fit, function(f) predict(f, type = type))
        expect_equal(unname(result$gradient), unname(numeric), tolerance = 2e-7)
        expect_true(all(is.finite(result$fit)))
      }
    }
})

test_that("population targets satisfy independence and rare-selection identities", {
  for (family in c("gaussian", "probit", "poisson", "pln")) {
    fit = prediction_fixture(family)
    fit$coefficients[grepl("^cor_", names(fit$coefficients))] = 0
    expect_equal(predict(fit, type = "selected_mean"), predict(fit), tolerance = 1e-12)
    expect_equal(predict(fit, type = "joint"),
                 predict(fit)*predict(fit, type = "selection_probability"), tolerance = 1e-12)
    fit$coefficients["selection.(Intercept)"] = -50
    rare = predict(fit, type = "selected_mean", se.fit = TRUE)
    # Independence keeps the selected mean equal to the response mean even
    # when selection is extremely rare.
    expect_equal(rare$fit, predict(fit), tolerance = 1e-10)
    expect_true(all(is.finite(rare$gradient)))
  }
})

test_that("history gradients include posterior weights and adaptive nodes", {
  for (family in c("gaussian", "probit", "poisson", "pln")) {
    fit = prediction_fixture(family, adaptive = TRUE)
    newdata = fit$data[1L, c("id", "x", "off"), drop = FALSE]
    history = fit$data[1:2, ]
    for (type in c("response", "selection_probability", "joint", "selected_mean")) {
      fun = function(f) predict(f, newdata, type = type,
        conditioning = "history", history = history)
      result = predict(fit, newdata, type = type, conditioning = "history",
                        history = history, se.fit = TRUE)
      expect_identical(colnames(result$gradient), names(fit$coefficients))
      expect_equal(unname(result$gradient),
                   unname(matrix(prediction_difference(fit, fun), 1)), tolerance = 3e-5)
    }
    reversed = predict(fit, newdata, conditioning = "history", history = history[2:1, ])
    expect_equal(reversed, predict(fit, newdata, conditioning = "history", history = history),
                 tolerance = 1e-10)
  }
})

test_that("history uses only explicit data and new subjects use the prior", {
  fit = prediction_fixture("probit")
  target = fit$data[1:2, ]
  history = fit$data[3:6, ]
  expect_equal(predict(fit, target, conditioning = "history", history = history),
               predict(fit, target), tolerance = 1e-12)
  target$y = 999; target$d = 999
  expect_equal(predict(fit, target, conditioning = "history", history = fit$data),
    predict(fit, fit$data[1:2, ], conditioning = "history", history = fit$data))
  expect_error(predict(fit, conditioning = "history"), "explicit history")
  expect_equal(predict(fit, conditioning = "history", history = fit$data[FALSE, ]), predict(fit))
})

test_that("prediction reconstructs factors, offsets, row order, and missing outcomes", {
  fit = prediction_fixture()
  data = fit$data
  data$f = factor(c("a", "b", "a", "b", "a", "b"))
  fit$prepared = prepare_data(y ~ f + offset(off), d ~ x, data, "id", "probit")
  names(fit$coefficients)[names(fit$coefficients) == "outcome.x"] = "outcome.fb"
  fit$data = data
  target = data[c(5, 1), c("id", "f", "x", "off")]
  expect_equal(predict(fit, target), predict(fit, data)[c(5, 1)])
  history = data[data$f == "a", ]; history$f = as.character(history$f)
  expect_true(all(is.finite(predict(fit, target, conditioning = "history", history = history))))
  bad = target; bad$f = "new"
  expect_error(predict(fit, bad), "new level")
  target$x[1] = NA_real_
  expect_true(is.na(predict(fit, target, type = "selection_probability")[1]))
  expect_true(all(is.finite(predict(fit, target))))
  expect_equal(predict(fit, target[, c("f", "off")]), predict(fit, target))
  expect_equal(predict(fit, target[, "x", drop=FALSE], type="selection_probability"),
    predict(fit, target, type="selection_probability"))
})

test_that("aggregate uncertainty averages gradients and both equation links include offsets", {
  fit = prediction_fixture()
  weights = c(1, 2, 3, 4, 5, 6)
  rows = predict(fit, se.fit = TRUE)
  average = predict(fit, se.fit = TRUE, aggregate = TRUE, weights = weights)
  expected = colSums(rows$gradient*weights/sum(weights))
  expect_equal(unname(average$gradient), matrix(unname(expected), 1))
  expect_equal(average$se.fit, sqrt(drop(expected %*% vcov(fit) %*% expected)))
  expect_equal(unname(predict(fit, type = "link")), .2 + .2*fit$data$x + fit$data$off)
  expect_equal(unname(predict(fit, type = "link", equation = "selection")), .2 + .2*fit$data$x)
  expect_error(predict(prediction_fixture(selection = FALSE), type = "link", equation = "selection"),
               "no selection")
  expect_error(predict(fit, aggregate = TRUE, weights = rep(0, 6)), "positive sum")
})

test_that("random effects retain outcome information from selection-only history", {
  fit = prediction_fixture("probit")
  fit$coefficients["cor_re"] = .5
  history = fit$data[fit$data$d == 0, ]
  effects = random_effects(fit, history)
  expect_identical(names(effects), c("id", "selection_mean", "selection_sd", "outcome_mean", "outcome_sd"))
  # Direct tensor integration with the same six-node rule; the outcome offset
  # enters the reference design with known coefficient one.
  theta = fit$coefficients
  par = c(alpha1 = theta[["selection.(Intercept)"]], alpha2 = theta[["selection.x"]],
    beta1 = theta[["outcome.(Intercept)"]], beta2 = theta[["outcome.x"]], beta3 = 1,
    theta[c("sd_re_selection", "sd_re_outcome", "cor_re", "cor_error")])
  for (k in seq_len(nrow(effects))) {
    rows = history[history$id == effects$id[k], ]
    panel = list(x = cbind(1, rows$x, rows$off), w = cbind(1, rows$x), y = rows$y, d = rows$d)
    expect_equal(unname(unlist(effects[k, -1L])),
      unname(ref_posterior_effects("selection_probit", par, panel, outer = 6L)), tolerance = 1e-10)
  }
  expect_true(any(abs(effects$outcome_mean) > 1e-4))
  expect_true(all(effects$outcome_sd > 0))
})

test_that("zero-column equations retain their offsets and latent effects", {
  data = data.frame(id=rep(1:3,each=2), y=c(0,1,1,0,0,1),d=c(1,1,1,0,1,0),off=seq(.1,.6,.1))
  for (selection in c(FALSE,TRUE)) {
    fit = prediction_fixture("probit",selection)
    fit$prepared = prepare_data(y~0+offset(off),if(selection) d~0 else NULL,data,"id","probit")
    fit$data = data
    fit$parameter_info = parameter_layout(fit$prepared)
    fit$coefficients = fit$coefficients[fit$parameter_info$name]
    fit$covariance = lapply(fit$covariance,function(V) V[names(fit$coefficients),names(fit$coefficients),drop=FALSE])
    expect_equal(unname(predict(fit)),pnorm(data$off/sqrt(1+.6^2)))
    expect_equal(unname(predict(fit,type="link")),data$off)
    if (selection) {
      expect_equal(unname(predict(fit,type="selection_probability")),rep(.5,6))
      expect_true(all(is.finite(predict(fit,data[1,,drop=FALSE],conditioning="history",history=data))))
    }
  }
})
