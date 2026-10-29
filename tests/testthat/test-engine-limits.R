# Engine work with no numerical effect: dimension-aware rule caches, lazy
# information steps and once-per-evaluator order checks; staged limit messages
# that name their controlling options; split GL capacity options and named byte
# bounds; PLN-only observation tolerances; and digest-bound staged meshes.

limits_data = function(seed = 3301, panels = 6L, waves = 2L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n))
  u = stats::rnorm(panels)
  data$d = as.integer(.4 + .5*data$x + .6*u[data$id] + stats::rnorm(n) > 0)
  data$y = stats::rpois(n, exp(.3 + .4*data$x + .5*u[data$id]))
  data$y[data$d == 0] = NA
  data
}

limits_theta = c("selection.(Intercept)" = .3, selection.x = .4, "outcome.(Intercept)" = .2,
  outcome.x = .3, sd_re_selection = .6, sd_re_outcome = .5, cor_re = .3)

test_that("outer rules are built only for dimensions a fit can use", {
  control = panel_control(nodes = 5L)
  zero = prepare_rule_control(control, 0L)$.rules
  full = prepare_rule_control(control)$.rules
  expect_length(zero$outer, 1L)
  expect_length(full$outer, 3L)
  expect_identical(zero$outer[[1L]], full$outer[[1L]])
  # A cached zero-dimension cache is extended when a panel rule is needed.
  upgraded = prepare_rule_control(prepare_rule_control(control, 0L), 2L)$.rules
  expect_identical(upgraded$outer, full$outer)
  expect_identical(prepare_rule_control(prepare_rule_control(control), 0L)$.rules, full)
})

test_that("information steps are computed only when observed information is", {
  data = limits_data()
  calls = 0L; original = information_steps
  local_mocked_bindings(information_steps = function(...) { calls <<- calls + 1L; original(...) },
    .package = "PanelSelect")
  fixed = limits_theta[c("sd_re_selection", "sd_re_outcome", "cor_re")]
  opg = panel_select(y ~ x, d ~ x, data, "id", "poisson", fixed = fixed,
    control = panel_control(nodes = 4L), init = "pooled")
  expect_identical(calls, 0L)
  observed = panel_select(y ~ x, d ~ x, data, "id", "poisson", fixed = fixed,
    control = panel_control(nodes = 4L), init = "pooled", vcov = "observed")
  expect_identical(calls, 1L)
  # A fresh evaluator starts adaptive modes cold, so agreement is to roundoff.
  expect_equal(vcov(opg, "observed"), vcov(observed, "observed"), tolerance = 1e-10)
})

test_that("an evaluator checks its parameter order once", {
  data = limits_data()
  prepared = prepare_data(y ~ x, d ~ x, data, "id", "poisson")
  info = parameter_layout(prepared)
  calls = 0L; original = check_parameter_order
  local_mocked_bindings(check_parameter_order = function(...) { calls <<- calls + 1L; original(...) },
    .package = "PanelSelect")
  evaluator = fit_evaluator(limits_theta, info, prepared, panel_control(nodes = 4L), numeric())
  expect_identical(calls, 1L)
  eta = parameter_to_optimizer(limits_theta, info)
  calls = 0L
  for (k in 1:3) evaluator$evaluate(eta + k*1e-3, trial = FALSE)
  expect_identical(calls, 0L)
  expect_error(fit_evaluator(rev(limits_theta), info, prepared, panel_control(), numeric()),
    "parameter order")
})

test_that("staged limit messages name the controlling option", {
  expect_identical(staged_options_name(panel_control(quadrature = "gl")), "quadrature_options")
  # An outer GH rule with a staged observation mesh uses the observation options.
  expect_identical(staged_options_name(panel_control(observation_quadrature = "gl")), "observation_options")
  expect_identical(staged_options_name(panel_control(observation_quadrature = "gl"), TRUE),
    "observation_options")
  # GL keeps its bound: the option to raise is the bound itself.
  expect_identical(staged_options_name(panel_control(quadrature = "gl"), bound = TRUE),
    "quadrature_options$bound")
  expect_identical(staged_options_name(panel_control(observation_quadrature = "gl"), bound = TRUE),
    "observation_options$bound")
  failure = tryCatch(staged_bound_failure("Outer omitted-tail envelope unresolved at maximum bound", 12,
    "quadrature_options$bound", c(a = 1), list(bound = 12)), panelselect_mesh_required = identity)
  expect_identical(conditionMessage(failure), paste("Staged quadrature: Outer omitted-tail envelope",
    "unresolved at maximum bound (maximum bound 12; raise quadrature_options$bound)"))
  expect_identical(failure$diagnostics, list(bound = 12))
  data = limits_data()
  budget = tryCatch(panel_select(y ~ x, d ~ x, data, "id", "poisson", fixed = limits_theta,
    control = panel_control(quadrature = "gl", quadrature_options = list(max_cells = 1))), error = identity)
  expect_s3_class(budget, "panelselect_integration_error")
  expect_match(conditionMessage(budget),
    "node or cell budget exhausted: [0-9]+ cells \\(quadrature_options\\$max_cells = 1\\)")
  bound = tryCatch(panel_select(y ~ x, d ~ x, data, "id", "poisson", fixed = limits_theta,
    control = panel_control(quadrature = "gl", quadrature_options = list(bound = 2))), error = identity)
  expect_match(conditionMessage(bound),
    "maximum bound \\(maximum bound 2; raise quadrature_options\\$bound\\)")
})

test_that("GL capacity options are explicit with the former defaults", {
  outer = panel_control(quadrature = "gl")
  expect_identical(outer$quadrature_options$max_cells, 4096L)
  expect_identical(outer$quadrature_options$batch_nodes, 8192L)
  inner = panel_control(observation_quadrature = "gl")
  expect_identical(inner$observation_options$batch_nodes, 128L)
  expect_error(panel_control(quadrature = "gl", quadrature_options = list(max_cells = 1.5)), "integer")
  set = panel_control(quadrature = "gl", quadrature_options = list(max_cells = 7, batch_nodes = 64))
  expect_identical(staged_options(set)[c("max_cells", "batch_nodes")],
    list(max_cells = 7L, batch_nodes = 64L))
  # Controls stored before these options existed resolve to the same limits.
  old = outer; old$quadrature_options[c("max_cells", "batch_nodes")] = NULL
  expect_identical(staged_options(old)[c("max_cells", "batch_nodes", "bound")],
    list(max_cells = 4096L, batch_nodes = 8192L, bound = 12))
  expect_identical(buffer_bytes, 16*1024^2)
  expect_identical(staged_buffer_bytes, 64*1024^2)
  expect_identical(staged_refinement_passes, 128L)
})

test_that("observation-rule tolerances gate information only for PLN", {
  control = panel_control(quadrature = "gl", observation_quadrature = "gl",
    observation_options = list(information_tol = 1e-6))
  expect_identical(staged_information_tolerance(control, "poisson"), 2e-3)
  expect_identical(staged_information_tolerance(control, "pln"), 1e-6)
})

test_that("a staged mesh is bound by a digest of its data arrays", {
  data = limits_data()
  fit = panel_select(y ~ x, d ~ x, data, "id", "poisson", fixed = limits_theta,
    control = panel_control(quadrature = "gl"))
  mesh = fit$control$integration_mesh
  expect_null(mesh$binding$x)
  expect_match(mesh$binding$data_md5, "^[0-9a-f]{32}$")
  prepared = fit$prepared; control = fit$control; fixed = fit$fixed
  expect_identical(staged_binding(prepared, control, fixed)$data_md5, mesh$binding$data_md5)
  for (field in c("x", "w", "y", "d", "offset_outcome", "offset_selection")) {
    changed = prepared
    changed[[field]][1L] = changed[[field]][1L] + 1
    expect_false(staged_binding_matches(mesh$binding, changed, control, fixed), label = field)
  }
  changed = prepared; changed$group = rev(changed$group)
  expect_false(staged_binding_matches(mesh$binding, changed, control, fixed))
  expect_true(staged_binding_matches(mesh$binding, unserialize(serialize(prepared, NULL)), control, fixed))
})
