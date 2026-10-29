# Output interface: summary and print lines, the exclusion-restriction field,
# strict method arguments, supplied covariance matrices, confint(parm =) forms,
# criteria sample checks, slimmer fitted objects and staged evaluation counts.

outputs_data = function(seed = 3101, panels = 16L, waves = 3L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    z = stats::rnorm(n), o = stats::runif(n, -.2, .2),
    region = factor(rep(c("a", "b", "c", "d"), length.out = n)))
  u = stats::rnorm(panels)
  data$d = as.integer(.3 + .5*data$x + .7*data$z + .6*u[data$id] + stats::rnorm(n) > 0)
  data$y = stats::rpois(n, exp(.2 + .4*data$x + .5*u[data$id]))
  data$y[data$d == 0] = NA
  data
}

# Complete starts with maxit = 0 give a cheap fitted object at known values.
outputs_fit = function(data, outcome, selection, family = "poisson", fixed = NULL,
                       id = "id", ...) {
  prepared = prepare_data(outcome, selection, data, id, family)
  info = parameter_layout(prepared)
  start = stats::setNames(ifelse(info$type == "sd", .5, ifelse(info$type == "correlation", .2, .1)),
    info$name)
  start = start[!names(start) %in% names(resolve_fixed(fixed, info, prepared))]
  panel_select(outcome, selection, data, id, family, start = start, fixed = fixed,
    control = panel_control(maxit = 0L, ...))
}

# The exclusion rule needs only the prepared designs and the fixed values.
exclusion_of = function(data, outcome, selection, fixed = NULL) {
  prepared = prepare_data(outcome, selection, data, "id", "poisson")
  info = parameter_layout(prepared)
  exclusion_restriction(list(prepared = prepared, fixed = resolve_fixed(fixed, info, prepared),
    data = data))
}

output_lines = function(x) capture.output(print(x))

test_that("summary and print name the model, sample, rule actually used and status", {
  data = outputs_data()
  both = outputs_fit(data, y ~ x, d ~ x + z, nodes = 5L)
  text = output_lines(summary(both))
  expect_identical(text[1L], "PanelSelect fit: Probit RE -> Poisson RE (panel selection)")
  expect_true("Outcome:   y ~ x" %in% text)
  expect_true("Selection: d ~ x + z" %in% text)
  expect_true("Sample:    48 rows in 16 panels (id \"id\")" %in% text)
  expect_true("Exclusion restriction: z" %in% text)
  expect_true("Integration: Fixed Gauss-Hermite (gh), 5 x 5 nodes (selection x outcome)" %in% text)
  expect_true(any(startsWith(text, "Optimizer:   not run (maxit = 0); Q_OPG = ")))
  expect_true(any(startsWith(text, "Covariance:  OPG (inverse sum of panel-score outer products), available")))
  expect_true(any(startsWith(text, "Correlation z and p values use the Fisher (atanh) scale; SD rows")))
  # One remaining effect takes its own element of a length-two order.
  outcome_only = outputs_fit(data, y ~ x, d ~ x + z, fixed = c(sd_re_selection = 0),
    quadrature = "gh", nodes = c(3L, 7L))
  text = output_lines(outcome_only)
  expect_true("Integration: Fixed Gauss-Hermite (gh), 7 nodes (outcome effect)" %in% text)
  expect_true("Fixed: sd_re_selection = 0; cor_re = 0 (required by sd_re_selection = 0)" %in% text)
  rows = output_lines(summary(outcome_only))
  expect_false(any(startsWith(rows, "sd_re_selection") | startsWith(rows, "cor_re")))
  # Cross-sectional selection needs no integral; the stored label says so.
  cross = outputs_fit(data, y ~ x, d ~ x + z, id = NULL, family = "linear")
  expect_identical(cross$integration$method, "none")
  text = output_lines(summary(cross))
  expect_identical(text[1L], "PanelSelect fit: Probit -> Linear (cross-sectional selection)")
  expect_true("Sample:    48 independent observations" %in% text)
  expect_true("Integration: None needed (closed-form likelihood)" %in% text)
  expect_true(any(grepl("BIC (observations)", text, fixed = TRUE)))
  # PLN reports its observation rule; a standalone model has no selection lines.
  pln = outputs_fit(data, y ~ x, d ~ x + z, "pln", nodes = 4L, observation_nodes = 6L)
  expect_true("Integration: Fixed Gauss-Hermite (gh), 4 x 4 nodes (selection x outcome); PLN observation rule: fixed Gauss-Hermite (gh), 6 nodes" %in%
    output_lines(pln))
  data$count = stats::rpois(nrow(data), 1)
  alone = outputs_fit(data, count ~ x, NULL, nodes = c(3L, 6L))
  text = output_lines(summary(alone))
  expect_identical(text[1L], "PanelSelect fit: Poisson RE (standalone panel)")
  expect_false(any(startsWith(text, "Selection:") | startsWith(text, "Exclusion")))
  expect_true("Integration: Fixed Gauss-Hermite (gh), 6 nodes (outcome effect)" %in% text)
  expect_null(summary(alone)$exclusion_restriction)
})

test_that("staged GL fits report the rule and mesh stages", {
  data = data.frame(id = rep(1:4, each = 2), d = c(1, 1, 1, 0, 0, 1, 1, 1),
    y = c(1, 0, 2, NA, NA, 1, 0, 3), x = c(-.5, .2, .4, -.1, .3, -.4, .1, .6))
  theta = c("selection.(Intercept)" = .2, "outcome.(Intercept)" = .1, outcome.x = .3,
    sd_re_selection = .5, sd_re_outcome = .4, cor_re = .2)
  fit = panel_select(y ~ x, d ~ 1, data, "id", "poisson", fixed = theta,
    control = panel_control(quadrature = "gl"))
  text = output_lines(fit)
  expect_true("Integration: Staged Gauss-Legendre (gl), order 20; 1 mesh stage used of at most 3" %in% text)
  expect_true("Optimizer:   not needed (all parameters fixed)" %in% text)
  expect_true("Exclusion restriction: none (every selection variable also enters the outcome equation)" %in%
    text)
})

test_that("the exclusion restriction uses raw variables, free columns and offsets", {
  data = outputs_data()
  expect_identical(exclusion_of(data, y ~ x, d ~ x + z), "z")
  expect_identical(exclusion_of(data, y ~ x, d ~ x), character())
  expect_identical(exclusion_of(data, y ~ x, d ~ poly(x, 2)), character())
  expect_identical(exclusion_of(data, y ~ x + z, d ~ x + z, c(outcome.z = 0)), "z")
  expect_identical(exclusion_of(data, y ~ x, d ~ x + z, c(selection.z = 0)), character())
  expect_identical(exclusion_of(data, y ~ x, d ~ x + offset(o)), "o")
  expect_identical(exclusion_of(data, y ~ x + z, d ~ x*z), character())
  expect_identical(exclusion_of(data, y ~ x, d ~ x:z), "z")
  expect_identical(exclusion_of(data, y ~ x, d ~ x + region, c(selection.regionb = 0)), "region")
  expect_identical(exclusion_of(data, y ~ x, d ~ x + region,
    c(selection.regionb = 0, selection.regionc = 0, selection.regiond = 0)), character())
  expect_identical(exclusion_of(data, y ~ x + o, d ~ .outcome + z), "z")
  # Fits saved before designs stored `assign` give the same answer.
  fit = outputs_fit(data, y ~ x, d ~ x + region + z, nodes = 4L)
  old = fit
  old$prepared$design$outcome$assign = old$prepared$design$selection$assign = NULL
  expect_identical(exclusion_restriction(old), c("region", "z"))
  expect_identical(exclusion_restriction(old), exclusion_restriction(fit))
  expect_identical(summary(fit)$exclusion_restriction, c("region", "z"))
})

test_that("methods reject misspelled arguments and accept covariance matrices", {
  data = outputs_data()
  fit = outputs_fit(data, y ~ x, d ~ x + z, nodes = 5L)
  expect_error(summary(fit, type = "observed"), "summary\\(\\) does not use the argument\\(s\\) type")
  expect_error(vcov(fit, se = "opg"), "selected with type = or vcov =")
  expect_error(confint(fit, type = "opg"), "confint\\(\\) does not use")
  expect_error(predict(fit, foo = 1), "predict\\(\\) does not use the argument\\(s\\) foo")
  V = vcov(fit, "opg")
  expect_equal(confint(fit, vcov = V), confint(fit, vcov = "opg"), tolerance = 1e-12)
  table = summary(fit, vcov = V)$coefficients
  free = !fit$parameter_info$fixed
  expect_equal(unname(table[free, "se"]), unname(sqrt(diag(V))[free]), tolerance = 1e-12)
  expect_identical(summary(fit, vcov = V)$covariance, "supplied")
  asymmetric = V; asymmetric[1L, 2L] = asymmetric[1L, 2L] + 1
  expect_error(confint(fit, vcov = asymmetric), "symmetric")
})

test_that("confint accepts negative and logical parameter selections", {
  data = outputs_data()
  fit = outputs_fit(data, y ~ x, d ~ x + z, nodes = 5L)
  p = length(coef(fit))
  expect_identical(confint(fit, parm = -1), confint(fit, parm = 2:p))
  expect_identical(confint(fit, parm = seq_len(p) %% 2 == 1), confint(fit, parm = seq(1L, p, 2L)))
  expect_identical(confint(fit, parm = "outcome.x"), confint(fit, parm = 5L))
  expect_error(confint(fit, parm = c(-1, 2)), "positive or negative indices")
  expect_error(confint(fit, parm = 0), "positive or negative indices")
  expect_error(confint(fit, parm = c(TRUE, FALSE)), "one TRUE or FALSE value per parameter")
  expect_error(confint(fit, parm = "nope"), "Unknown parameter in parm: nope")
})

test_that("criteria warn when compared fits use different samples", {
  data = outputs_data()
  fit = outputs_fit(data, y ~ x, d ~ x + z, nodes = 4L)
  same = outputs_fit(data, y ~ x, d ~ x, nodes = 4L)
  fewer = outputs_fit(data[data$id <= 12L, ], y ~ x, d ~ x + z, nodes = 4L)
  expect_silent(AIC(fit, same))
  expect_silent(BIC(fit, same))
  expect_warning(AIC(fit, fewer), "not all fitted to the same sample")
  expect_warning(BIC(fit, fewer, unit = "observations"), "not all fitted to the same sample")
  expect_silent(AIC(fewer))
})

test_that("residual types and slimmer fitted objects are explicit", {
  data = outputs_data()
  fit = outputs_fit(data, y ~ x, d ~ x + z, nodes = 4L)
  expect_error(residuals(fit, type = "pearson"), '"selected_mean".*"pearson"')
  expect_identical(residuals(fit, type = "resp"), residuals(fit, type = "response"))
  expect_null(fit$sample)
  expect_null(fit$optimizer$hessian)
  expect_null(fit$optimizer$rejected)
  expect_identical(rownames(fit$score), as.character(1:16))
  expect_identical(names(fit$loglik_by_unit), as.character(1:16))
  expect_identical(sum(fit$loglik_by_unit), fit$loglik)
})

test_that("the boundary note names free parameters near their boundary", {
  model = resolve_model("gaussian", TRUE, "panel")
  info = data.frame(name = c("outcome.x", "sd_re_selection", "sd_re_outcome", "sd_error", "cor_re", "cor_error"),
    type = c("coefficient", "sd", "sd", "sd", "correlation", "correlation"), fixed = FALSE)
  theta = c(outcome.x = .5, sd_re_selection = .8, sd_re_outcome = .3, sd_error = 2, cor_re = .4, cor_error = .2)
  expect_null(boundary_note(theta, info, model))
  near = theta; near[c("cor_re", "sd_re_selection")] = c(.9995, 4e-4)
  expect_identical(boundary_note(near, info, model), paste("sd_re_selection = 4e-04, cor_re = 0.9995 are",
    "near the parameter boundary, where this screen can overstate the remaining distance"))
  # A Linear outcome SD is judged against the residual SD; fixed values are not free.
  scaled = theta; scaled[c("sd_re_outcome", "sd_error")] = c(.01, 1000)
  expect_match(boundary_note(scaled, info, model), "^sd_re_outcome = 0.01 is near")
  fixed = info; fixed$fixed[fixed$name == "cor_re"] = TRUE
  at_fixed = theta; at_fixed["cor_re"] = .9995
  expect_null(boundary_note(at_fixed, fixed, model))
})

test_that("fitted contrasts survive re-framing without warnings", {
  data = outputs_data()
  contrasts(data$region) = stats::contr.sum(4)
  fit = outputs_fit(data, y ~ x + region, d ~ x + z, nodes = 4L)
  expect_silent(matrix <- model.matrix(fit))
  expected = stats::model.matrix(~ x + region, data)
  expect_identical(colnames(matrix), colnames(expected))
  expect_equal(unname(matrix), unname(expected))
  expect_identical(attr(matrix, "contrasts"), attr(expected, "contrasts"))
  expect_silent(predict(fit, newdata = data[1:3, ]))
})

test_that("old fitted objects and saved summaries still print", {
  data = outputs_data()
  fit = outputs_fit(data, y ~ x, d ~ x + z, nodes = 4L)
  old = fit
  old$prepared$design$outcome$assign = old$prepared$design$selection$assign = NULL
  old$optimizer$evaluations = NULL
  old$sample = list(n_observations = 48L, n_units = 16L, n_panels = 16L)
  old$optimizer$hessian = -old$optimizer$information
  names(old$loglik_by_unit) = NULL; rownames(old$score) = NULL
  expect_identical(output_lines(summary(old)), output_lines(summary(fit)))
  saved = structure(list(model_id = "probit_re_poisson_re", structure = "panel",
    coefficients = summary(fit)$coefficients, fixed = character(), covariance = "opg",
    n_rows = 48L, n_panels = 16L, loglik = fit$loglik, AIC = AIC(fit), BIC = BIC(fit),
    converged = FALSE), class = "summary.panelselect")
  text = output_lines(saved)
  expect_identical(text[1L], "PanelSelect probit_re_poisson_re model")
  expect_true(any(startsWith(text, "Covariance:  OPG (inverse sum of panel-score outer products)")))
})

test_that("staged evaluation counts cover every stage", {
  data = outputs_data()
  fit = panel_select(y ~ x, d ~ x + z, data, "id", "poisson", fixed = c(cor_re = .3),
    control = panel_control(nodes = 4L), init = "pooled")
  expect_gte(fit$optimizer$evaluations, fit$optimizer$counts[["function"]])
  staged = fit
  staged$integration$staging = list(stages_used = 2L, maximum_stages = 3L)
  staged$optimizer$evaluations = 30L
  expect_match(optimizer_description(staged), "in the final mesh stage; 30 likelihood evaluations over 2 stages")
})
