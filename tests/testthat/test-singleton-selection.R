singleton_selection_data = function() data.frame(id = 1:8,
  x = c(-1, -.6, -.2, .1, .4, .7, 1, 1.3),
  z = c(0, 1, 0, -1, 1, 0, -1, 1),
  d = c(0, 1, 1, 0, 1, 0, 1, 1), y = c(NA, 0, 1, NA, 0, NA, 1, 0), off = 0)

test_that("unrestricted singleton selection models reject before initialization", {
  data = singleton_selection_data()
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  for (family in c("probit", "gaussian", "poisson", "pln")) {
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family),
      "unidentified singleton likelihood ridge")
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family,
      start = c(sd_re_selection = .7, sd_re_outcome = .6, cor_re = 0)),
      "unidentified singleton likelihood ridge")
  }
})

test_that("free correlations and two free angles produce the demonstrated ridges", {
  data = singleton_selection_data()
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  for (family in c("probit", "gaussian", "pln")) {
    fixed = c(sd_re_selection = .7, sd_re_outcome = .6,
      if (family != "probit") c(sd_error = .8))
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family,
      fixed = fixed), "cor_re, cor_error")
    # Both scale angles can move even at the equal-correlation critical diagonal.
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family,
      fixed = c(cor_re = .3, cor_error = .3)), "unidentified singleton")
  }
})

test_that("known indices and fixed parameters retain informative singleton exceptions", {
  data = singleton_selection_data(); data$off = data$z
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  # An outside-span selection index anchors its scale, leaving one free correlation.
  expect_error(panel_select(y ~ x, d ~ x + offset(off), data, "id", "poisson"),
    "Initialization reached")
  expect_error(panel_select(y ~ x, d ~ x + z, data, "id", "poisson",
    fixed = c(selection.z = 1)), "Initialization reached")
  # Fixed nonzero person correlation identifies the sole available selection angle.
  expect_error(panel_select(y ~ x, d ~ x, data, "id", "poisson",
    fixed = c(cor_re = .3)), "Initialization reached")
  # A fixed nonzero correlation can distinguish an otherwise free variance split.
  for (family in c("gaussian", "pln")) {
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family,
      fixed = c(sd_re_selection = .7, cor_re = .3, cor_error = 0)), "Initialization reached")
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family,
      fixed = c(sd_re_selection = 0, cor_error = .3)), "Initialization reached")
  }
  # Both probit indices can be anchored while just one correlation remains free.
  expect_error(panel_select(y ~ x + offset(off),
    d ~ x + offset(off), data, "id", "probit", fixed = c(cor_error = 0)), "Initialization reached")
  for (family in c("probit", "gaussian", "poisson", "pln")) {
    prepared = prepare_data(y ~ x, d ~ x, data, "id", family)
    info = parameter_layout(prepared)
    fixed = setNames(rep(.2, nrow(info)), info$name); info$fixed = TRUE
    expect_silent(validate_design_rank(prepared, info, fixed))
  }
})

test_that("only exact fixed zeros remove correlation information", {
  data = singleton_selection_data()
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  expect_error(panel_select(y ~ x, d ~ x, data, "id", "poisson",
    fixed = c(cor_re = 0)), "unidentified singleton")
  expect_error(panel_select(y ~ x, d ~ x, data, "id", "poisson",
    fixed = c(cor_re = 1e-100)), "Initialization reached")
  for (family in c("gaussian", "pln"))
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family,
      fixed = c(sd_re_selection = 0, cor_error = 0)), "sd_re_outcome/sd_error split")
  expect_error(panel_select(y ~ x, d ~ x, data, "id", "pln",
    fixed = c(sd_re_outcome = 0, sd_error = 0)), "unidentified singleton")
  # Zero outcome variance is allowed once an independent known index anchors selection.
  data$off = data$z
  expect_error(panel_select(y ~ x, d ~ x + offset(off), data, "id", "pln",
    fixed = c(sd_re_outcome = 0, sd_error = 0)), "Initialization reached")
})

test_that("span checks use actual fixed values and observed outcome rows", {
  data = singleton_selection_data()
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  data$off = .3 * data$x
  expect_error(panel_select(y ~ x, d ~ x + offset(off), data, "id", "poisson"),
    "unidentified singleton")
  data$off = -data$z
  expect_error(panel_select(y ~ x, d ~ x + z + offset(off), data, "id", "poisson",
    fixed = c(selection.z = 1)), "unidentified singleton")
  # An offset only on unselected rows does not anchor the observed outcome scale.
  data$off = ifelse(data$d == 0, data$z, 0)
  expect_error(panel_select(y ~ x + offset(off), d ~ x, data, "id", "probit",
    fixed = c(sd_re_selection = .7, cor_re = 0, cor_error = 0)), "unidentified singleton")
  data$off = data$z
  expect_error(panel_select(y ~ x + offset(off), d ~ x, data, "id", "probit",
    fixed = c(sd_re_selection = .7, cor_re = 0, cor_error = 0)), "Initialization reached")
  data$huge = 2
  expect_error(panel_select(y ~ x, d ~ x + huge, data, "id", "poisson",
    fixed = c(selection.huge = 1e308)), "fixed selection index is not finite")
})

test_that("existing diagnostics run first and retained panel lengths determine scope", {
  data = singleton_selection_data()
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  data$duplicate = data$x
  expect_error(panel_select(y ~ x, d ~ x + duplicate, data, "id", "probit"),
    "free selection coefficients are not identified")
  data$d = 1; data$y = rep(0:1, 4)
  expect_error(panel_select(y ~ x, d ~ x, data, "id", "probit"),
    "no finite maximum-likelihood")
  data$d = 0; data$y = NA_real_
  expect_error(panel_select(y ~ 0, d ~ 0, data, "id", "probit"),
    "No outcomes are observed")
  data = singleton_selection_data(); data$id[2] = data$id[1]
  for (family in c("probit", "gaussian", "poisson", "pln"))
    expect_error(panel_select(y ~ x, d ~ x, data, "id", family), "Initialization reached")
  omitted = rbind(singleton_selection_data(), singleton_selection_data())
  omitted$d[9:16] = NA
  expect_error(panel_select(y ~ x, d ~ x, omitted, "id", "poisson", na_action = "omit"), "unidentified singleton")
})

test_that("starts do not replace standalone or joint identification restrictions", {
  data = singleton_selection_data()
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  expect_error(panel_select(y ~ x, d ~ x, data, "id", "probit", init = "zero",
    start = c(sd_re_selection = .7, sd_re_outcome = .6, cor_re = 0, cor_error = 0)),
    "unidentified singleton")
  expect_error(panel_select(d ~ x, NULL, data, "id", "probit", start = c(sd_re_outcome = .7)),
    "standalone probit")
})
