test_that("singleton PLN rejects two free scales before initialization", {
  data = data.frame(id = 1:6, x = c(-1, -.4, 0, .3, .7, 1.2),
                     y = c(0, 1, 3, 2, 4, 1), off = seq(-.2, .3, length.out = 6))
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  for (fixed in list(NULL, c("outcome.(Intercept)" = .2, outcome.x = -.3))) {
    expect_error(panel_select(y ~ x + offset(off), NULL, data, "id", "pln", fixed = fixed),
      "cannot identify both sd_re_outcome and sd_error; fix one scale or supply repeated")
  }
})

test_that("the singleton PLN guard uses retained rows after omission", {
  data = data.frame(id = rep(1:3, each = 2), x = rep(-1:1, each = 2),
                     y = c(1, NA, 2, NA, 0, NA))
  local_mocked_bindings(initial_parameters = function(...) stop("Initialization reached"))
  expect_error(panel_select(y ~ x, NULL, data, "id", "pln", na_action = "omit"), "one retained row per panel")
})

test_that("fixing either PLN scale preserves singleton input support", {
  data = data.frame(id = 1:6, x = c(-1, -.4, 0, .3, .7, 1.2), y = c(0, 1, 3, 2, 4, 1))
  prepared = prepare_data(y ~ x, NULL, data, "id", "pln")
  for (scale in c("sd_re_outcome", "sd_error")) for (value in c(0, .7)) {
    info = parameter_layout(prepared)
    fixed = resolve_fixed(setNames(value, scale), info, prepared)
    info$fixed = info$name %in% names(fixed)
    expect_silent(validate_design_rank(prepared, info))
  }
  info = parameter_layout(prepared)
  info$fixed = TRUE
  expect_silent(validate_design_rank(prepared, info))
})

test_that("repeated panels and nested singleton models use their respective guards", {
  data = data.frame(id = c(1, 1, 2, 3, 4, 5), x = c(-1, -.4, 0, .3, .7, 1.2),
                     y = c(0, 1, 3, 2, 4, 1), d = c(1, 0, 1, 0, 1, 0))
  repeated = prepare_data(y ~ x, NULL, data, "id", "pln")
  expect_silent(validate_design_rank(repeated, parameter_layout(repeated)))
  data$id = seq_len(nrow(data))
  selected = prepare_data(y ~ x, d ~ x, data, "id", "pln")
  expect_error(validate_design_rank(selected, parameter_layout(selected)), "unidentified singleton")
  poisson = prepare_data(y ~ x, NULL, data, "id", "poisson")
  expect_silent(validate_design_rank(poisson, parameter_layout(poisson)))
})
