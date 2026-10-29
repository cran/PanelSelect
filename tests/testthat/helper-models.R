# Independent fixtures bypass formula preparation so numerical failures are local.
#
# Bitwise checks: expect_identical() between two production paths states a
# contract of this package's arithmetic, not a tolerance. The paths perform the
# same floating-point operations in the same order, on one thread with the
# reference BLAS. Other compilers, FMA contraction or a threaded BLAS can break
# such a check without a defect. Comparisons with independent references always
# use explicit tolerances.
reference_path = testthat::test_path('fixtures', 'likelihoods.R')
stopifnot(file.exists(reference_path))

# Frozen fixtures keep their original parameter names (sd_selection and
# sd_outcome); names, dimnames and character values are
# translated when a fixture is loaded, never in the fixture file.
fixture_names = function(x) {
  if (is.list(x)) for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] = fixture_names(x[[i]])
  if (is.character(x)) x = current_names(x)
  if (!is.null(names(x))) names(x) = current_names(names(x))
  if (!is.null(dim(x)) && !is.null(dimnames(x)))
    dimnames(x) = lapply(dimnames(x), function(d) if (is.null(d)) d else current_names(d))
  x
}
core_references = fixture_names(dget(reference_path))

core_models = c('probit_re', 'poisson_re', 'pln_re', 'selection_gaussian',
                 'selection_probit', 'selection_poisson', 'selection_pln')

core_control = function(method = 'gh', nodes = 7L, inner = 13L, exact_reductions = TRUE) {
  list(quadrature = method, nodes = nodes, observation_nodes = inner,
       exact_reductions = exact_reductions, mode_maxit = 60L, mode_tol = 1e-10)
}

core_fixture = function(model, correlation_sign = 1) {
  selection = startsWith(model, 'selection_')
  family = sub('^selection_|_re$', '', model)
  # Explicit panel metadata, as prepare_data() records it; never inferred from absence.
  model_id = c(probit_re = 'probit_re', poisson_re = 'poisson_re', pln_re = 'pln_re',
    selection_gaussian = 'probit_re_linear_re', selection_probit = 'probit_re_probit_re',
    selection_poisson = 'probit_re_poisson_re', selection_pln = 'probit_re_pln_re')[[model]]
  x = cbind('(Intercept)' = 1, x = c(-0.4, 0.3, 0.7))
  w = cbind('(Intercept)' = 1, w = c(0.5, -0.3, 0.1))
  y = switch(family, probit = c(1, 0, 0), gaussian = c(0.2, -0.1, -0.4), c(1, 0, 2))
  d = if (selection) c(1, 0, 1) else rep(1, 3)
  if (selection) y[d == 0] = NA_real_
  theta = c('outcome.(Intercept)' = 0.25, outcome.x = -0.2, sd_re_outcome = 0.4)
  if (selection) theta = c('selection.(Intercept)' = -0.15, selection.w = 0.3,
    theta, sd_re_selection = 0.45, cor_re = -0.35*correlation_sign)
  if (family %in% c('gaussian', 'pln')) theta = c(theta, sd_error = 0.5)
  if (selection && family != 'poisson') theta = c(theta, cor_error = 0.25*correlation_sign)
  # Canonical parameter_layout() order, on which positional optimizer maps rely.
  theta = theta[c(grep('^selection\\.', names(theta), value = TRUE),
    grep('^outcome\\.', names(theta), value = TRUE),
    intersect(c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error'), names(theta)))]
  prepared = list(model = model, model_id = model_id, structure = 'panel',
    family = family, has_selection = selection,
    x = x, w = if (selection) w else matrix(numeric(), 3, 0), y = y, d = d,
    offset_outcome = rep(0, 3), offset_selection = rep(0, 3), group = rep(1L, 3),
    group_ptr = c(0L, 3L), id_values = 'person', n_rows = 3L, n_panels = 1L,
    selected = d == 1)
  stopifnot(identical(names(theta), parameter_layout(prepared)$name))
  reference = theta
  names(reference)[startsWith(names(reference), 'outcome.')] = c('beta1', 'beta2')
  if (selection) names(reference)[startsWith(names(reference), 'selection.')] = c('alpha1', 'alpha2')
  list(theta = theta, prepared = prepared, reference = reference,
       panel = list(x = x, w = w, y = y, d = d))
}

# The frozen likelihood fixture keeps its original hand-built parameter order.
# Index its theta and natural score by name to the fixture's canonical order.
core_reference = function(fixture, reference) {
  stopifnot(setequal(names(reference$theta), names(fixture$theta)))
  score = stats::setNames(reference$score, names(reference$theta))
  reference$theta = reference$theta[names(fixture$theta)]
  reference$score = unname(score[names(fixture$theta)])
  reference
}

core_difference = function(fun, theta, step = 1e-5) {
  vapply(seq_along(theta), function(j) {
    h = step*max(1, abs(theta[j]))
    plus = minus = theta; plus[j] = plus[j]+h; minus[j] = minus[j]-h
    (fun(plus)-fun(minus))/(2*h)
  }, numeric(1))
}

expect_core_score = function(fixture, control, fixed = numeric(), tolerance = 3e-6) {
  actual = evaluate_model(fixture$theta, fixture$prepared, control, fixed)
  numerical = core_difference(function(theta) sum(evaluate_model(theta,
    fixture$prepared, control, fixed, score = FALSE)$loglik), fixture$theta)
  free = !names(fixture$theta) %in% names(fixed)
  expect_lt(max(abs(colSums(actual$score)[free]-numerical[free])), tolerance)
  if (any(!free)) expect_equal(unname(actual$score[, !free]), rep(0, sum(!free)))
  invisible(actual)
}

# Value and warning messages of an expression, so a test can assert exactly the
# warnings a call is known to give instead of suppressing every warning.
with_warnings = function(expr) {
  messages = character()
  value = withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart('muffleWarning')
  })
  list(value = value, warnings = messages)
}
