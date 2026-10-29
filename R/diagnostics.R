# 1. Public integration sensitivity check

#' Check integration sensitivity at the fitted estimates
#'
#' Re-evaluate a fitted model using higher quadrature orders or an alternative GL
#' rule. The comparisons hold parameters fixed, so they isolate changes in the
#' numerical likelihood and score. A higher order is not guaranteed to be more
#' accurate, and agreement is evidence of stability rather than a certificate.
#'
#' @param object A fitted `panelselect` model. Its retained data, restrictions
#'   and parameters define every comparison; the input object is not modified.
#' @param nodes Person-effect orders for the comparison. With
#'   `compare = "higher_order"`, the default adds the larger of 5 or half the
#'   fitted order (rounded up) to each fitted order. A scalar uses the same order
#'   on both effects; a pair follows selection then outcome. Supplied orders
#'   cannot decrease. With `compare = "gl"`, this is the GL order per segment
#'   (default 20), not a total node budget.
#' @param observation_nodes PLN observation order for the comparison: twice the
#'   fitted order by default for `"higher_order"`, or GL order 20 for `"gl"`.
#'   Only active PLN observation integration contributes a comparison row.
#' @param tolerance Named positive numeric vector with `loglik` (default `1e-5`)
#'   and `score` (default `1e-4`). The log-likelihood tolerance applies to the
#'   absolute total change and the largest absolute independent-unit change;
#'   the score tolerance applies to the largest absolute total and unit
#'   optimizer-score component changes. These are comparison tolerances, not
#'   optimizer stopping controls or certified integral-error bounds.
#' @param refit Logical; default `FALSE` performs only fixed-parameter checks.
#'   `TRUE` additionally calls [update()] with the compared controls and fitted
#'   natural coefficients as warm starts, retaining the resulting fit in
#'   `$refit`. Fitting errors are not silently replaced by a different method.
#' @param compare `"higher_order"` (default) refines the fitted rule;
#'   `"gl"` compares a GH/AGH person-effect fit against a separately constructed
#'   GL rule, also using GL for an active PLN observation integral. Fits already
#'   using GL should use `"higher_order"`. See the comparison table below.
#' @param quadrature_options Named GL options, with names and meanings from
#'   [panel_control()], used only for `compare = "gl"`. Defaults are obtained
#'   from a GL control; increase bounds or budgets explicitly when a comparison
#'   reports the corresponding limit.
#' @param ... Unused by the print method. `check_quadrature()` rejects unknown
#'   arguments and gives replacements for earlier argument names.
#'
#' @section Choose a comparison:
#' | Option | What changes | When it helps |
#' |---|---|---|
#' | `"higher_order"` | Orders increase within the fitted integration methods. Staged rules can require refined meshes. | First check for order sensitivity; also the supported check of an existing GL fit. |
#' | `"gl"` | GH/AGH person-effect integration is compared with GL at the same estimates, using moving cuts or a fresh staged mesh as appropriate. | Investigate method sensitivity, especially large scales, rare histories, sharp transitions or unstable estimates. |
#'
#' For PLN fits, `person_effect` changes only the outer integral, `observation`
#' only the occasion-level integral, and `both` changes both where both are
#' active. Fixed restrictions and exact reductions can remove an integral, so
#' only effective changes produce rows. A fully closed-form model has no
#' numerical integral to refine and the function stops with that explanation.
#'
#' @section Read the results and decide what to do:
#' | Comparison column | Meaning |
#' |---|---|
#' | `rule` | `person_effect`, `observation`, `both`, or `gl`. |
#' | `nodes`, `reference_nodes` | Compared and fitted person-effect orders, respectively. |
#' | `observation_nodes`, `observation_quadrature` | Compared PLN inner rule when applicable. |
#' | `mesh_changed` | Whether a staged comparison changed the retained mesh; `NA` when not applicable. |
#' | `loglik_change` | Signed change in the total log likelihood. |
#' | `max_panel_change` | Largest absolute independent-unit log-likelihood change. |
#' | `max_score_change`, `max_panel_score_change` | Largest absolute total and individual-unit optimizer-score changes. |
#' | `integration_failed`, `message` | Numerical failure and its explanation, separate from tolerance exceedance. |
#' | `within_tolerance` | Whether all four discrepancies meet their tolerances; `NA` after numerical failure. |
#'
#' - If discrepancies exceed tolerance, investigate larger orders or another
#'   method and assess changes in estimates/uncertainty after an explicit refit.
#' - If integration fails, inspect `message`, `$failures` and `$integration`;
#'   change only the named controls whose limits or accuracy require attention.
#' - Even small changes do not prove both rules accurate. Parameter-dependent
#'   integration error can move an optimum despite correct finite-rule scores.
#'   Compare competing endpoints under a common checked rule when needed.
#' - Keep optimizer status, [check_stationarity()], covariance availability,
#'   and substantive sensitivity separate from this integration comparison.
#'
#' @return A `panelselect_quadrature_check` list with `compare`, `fitted` (a
#'   fitted-rule description), the `comparison` data frame, `tolerance`,
#'   `optimizer_converged`, `covariance_valid`, `integration`, `failures`,
#'   `refit` (`NULL` by default), and `qualification`. `print()` displays the
#'   main comparison columns and returns the object invisibly. Covariance
#'   validity describes the requested fitted estimator, not statistical calibration.
#' @family diagnostic functions
#' @seealso [panel_control()], [check_stationarity()], [update.panelselect()]
#' @examples
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#' fit = panel_select(y ~ x, observed ~ x + z, dat, "person", "linear")
#' check_quadrature(fit)
#' check_quadrature(fit, nodes = 20, tolerance = c(loglik = 1e-6, score = 1e-5))
#' @export
check_quadrature = function(object, nodes = NULL, observation_nodes = NULL,
                             tolerance = c(loglik = 1e-5, score = 1e-4),
                             refit = FALSE, compare = c("higher_order", "gl"),
                             quadrature_options = list(), ...) {
  dots = list(...)
  renamed = c(method = "compare", inner_nodes = "observation_nodes", gk = "quadrature_options")
  earlier = intersect(names(dots), names(renamed))
  if (length(earlier))
    stop(paste0(earlier, " is now ", renamed[earlier], collapse = "; "),
         if ("method" %in% earlier) '; method = "refine" is now compare = "higher_order"', call. = FALSE)
  reject_unused(dots, "check_quadrature", "see ?check_quadrature")
  if (!inherits(object, "panelselect")) stop("object must be a panelselect fit", call. = FALSE)
  object = current_fit(object)
  stop_removed_settings(object$control)
  # Legacy optimizer options matter only to a refit, which stops on them.
  if (!refit) attr(object$control, "legacy_optimizer_control") = NULL
  if (!is.numeric(tolerance) || length(tolerance) != 2L ||
      !setequal(names(tolerance), c("loglik", "score")) ||
      any(!is.finite(tolerance) | tolerance <= 0))
    stop("tolerance must name positive loglik and score tolerances", call. = FALSE)
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit))
    stop("refit must be TRUE or FALSE", call. = FALSE)
  if (identical(compare, "refine")) stop('compare = "refine" is now "higher_order"', call. = FALSE)
  if (identical(compare, "gk")) stop('compare = "gk" was removed; use compare = "gl"', call. = FALSE)
  compare = match.arg(compare)
  control = object$control
  if (compare == "gl")
    return(check_gl_reference(object, nodes, observation_nodes, tolerance, refit, quadrature_options))
  if (length(quadrature_options)) stop("quadrature_options apply to compare = 'gl'", call. = FALSE)
  if (uses_staged_integration(object$prepared, control, object$fixed))
    return(check_staged_quadrature(object, nodes, observation_nodes, tolerance, refit))
  if (uses_gl_quadrature(control))
    return(check_split_quadrature(object, nodes, observation_nodes, tolerance, refit))
  nodes = nodes %||% (control$nodes + pmax(5L, ceiling(control$nodes/2)))
  observation_nodes = observation_nodes %||% (2L*control$observation_nodes)
  refined = control; refined$nodes = nodes; refined$observation_nodes = observation_nodes
  refined = validate_panel_control(refined)
  if (any(rep_len(refined$nodes, 2L) < rep_len(control$nodes, 2L)) ||
      refined$observation_nodes < control$observation_nodes)
    stop("Refined orders must not be smaller than fitted orders", call. = FALSE)
  theta = coef(object)
  baseline = evaluate_model(theta, object$prepared, control, object$fixed)
  outer_orders = quadrature_outer_orders(baseline$diagnostics)
  inner_active = object$prepared$family == "pln" && theta[["sd_error"]] > 0 &&
    any(object$prepared$d == 1)
  free = !object$parameter_info$fixed
  jacobian = parameter_jacobian(theta, object$parameter_info)[free]
  # Rows: person_effect and observation refine one rule each, both refines both.
  variants = list()
  if (length(outer_orders)) {
    variants$person_effect = control
    variants$person_effect$nodes = refined$nodes
  }
  if (inner_active) {
    variants$observation = control; variants$observation$observation_nodes = refined$observation_nodes
  }
  changed = function(current) !identical(rep_len(current$nodes, 2L), rep_len(control$nodes, 2L)) ||
    (inner_active && current$observation_nodes != control$observation_nodes)
  variants = Filter(changed, variants)
  if (!length(variants)) {
    if (!length(outer_orders) && !inner_active)
      stop("This model has no numerical integral to refine", call. = FALSE)
    stop("At least one effective quadrature order must increase", call. = FALSE)
  }
  compare_variant = function(name) {
    current = variants[[name]]
    result = tryCatch(evaluate_model(theta, object$prepared, current, object$fixed),
                       error = identity)
    failed = inherits(result, "error")
    # Dimension reductions and independent factors may use only one supplied
    # axis. Compare the rules actually evaluated, not just raw control vectors.
    if (name == "person_effect" && !failed &&
        identical(quadrature_outer_orders(result$diagnostics), outer_orders))
      return(NULL)
    comparison_row(name, current$nodes, control$nodes,
      if (failed) rep(NA_real_, 4L) else refinement_change(result, baseline, free, jacobian),
      failed, if (failed) conditionMessage(result) else "", tolerance,
      observation_nodes = if (inner_active) current$observation_nodes else NA_integer_,
      observation_quadrature = if (inner_active) current$observation_quadrature %||% "gh" else NA_character_)
  }
  comparison = do.call(rbind, lapply(names(variants), compare_variant))
  if (is.null(comparison)) stop("At least one effective quadrature order must increase", call. = FALSE)
  if (all(c("person_effect", "observation") %in% comparison$rule)) {
    variants$both = refined
    comparison = rbind(comparison, compare_variant("both"))
  }
  # Refits go through update(): warm starts, subset rows and the call record.
  fitted = if (refit) do.call(stats::update, list(object, control = refined)) else NULL
  quadrature_check(object, "higher_order", comparison, tolerance, fitted,
    qualification = "Same-method order refinement at the fitted parameters; no accuracy certificate or silent refit.")
}

# 2. Staged, split and independent-rule comparisons

check_staged_quadrature = function(object, nodes, observation_nodes, tolerance, refit) {
  control = object$control; theta = coef(object)
  baseline = evaluate_model(theta, object$prepared, control, object$fixed)
  native_outer = control$quadrature %in% c('gh', 'agh') && !uses_gl_quadrature(control)
  outer_orders = if (native_outer) quadrature_outer_orders(baseline$diagnostics$native) else integer()
  outer = if (native_outer) length(outer_orders) > 0L else
    any(vapply(baseline$diagnostics$panels, function(x) x$dimension > 0L, logical(1L)))
  inner = object$prepared$family == 'pln' && theta[['sd_error']] > 0 &&
    any(object$prepared$d == 1 & (object$prepared$active %||% TRUE))
  if (!outer && !inner) stop('This model has no numerical integral to refine', call. = FALSE)
  refined = control; variants = list()
  if (outer) {
    refined$nodes = nodes %||% (control$nodes + pmax(5L, ceiling(control$nodes/2)))
    # Compare values: a supplied order may be double where the fitted one is integer.
    orders = rep_len(refined$nodes, 2L); fitted_orders = rep_len(control$nodes, 2L)
    if (any(orders < fitted_orders) || all(orders == fitted_orders))
      stop('At least one effective person-effect order must increase and none may decrease', call. = FALSE)
    if (identical(control$quadrature, 'gl')) refined$quadrature_options$reference_nodes = control$nodes
    variants$person_effect = control
    variants$person_effect$nodes = refined$nodes
    variants$person_effect$quadrature_options = refined$quadrature_options
  } else if (!is.null(nodes)) stop('This model has no person-effect integral to refine', call. = FALSE)
  if (inner) {
    refined$observation_nodes = observation_nodes %||% (2L*control$observation_nodes)
    if (refined$observation_nodes <= control$observation_nodes)
      stop('The effective observation order must increase', call. = FALSE)
    if (identical(control$observation_quadrature, 'gl'))
      refined$observation_options$reference_nodes = control$observation_nodes
    variants$observation = control; variants$observation$observation_nodes = refined$observation_nodes
    variants$observation$observation_options = refined$observation_options
  } else if (!is.null(observation_nodes))
    stop('This model has no PLN observation integral to refine', call. = FALSE)
  if (outer && inner) variants$both = refined
  free = !object$parameter_info$fixed
  jacobian = parameter_jacobian(theta, object$parameter_info)[free]
  fitted_mesh = control$integration_mesh
  results = failures = list()
  comparison = lapply(names(variants), function(name) {
    current = variants[[name]]; current$integration_mesh = current$.rules = NULL
    # An invalid requested order or option stops; only integration may fail here.
    current = validate_panel_control(current)
    result = tryCatch({
      current$integration_mesh = prepare_staged_integration(theta, object$prepared,
        current, object$fixed, previous = fitted_mesh, epoch = 1L)
      validate_staged_integration(theta, object$prepared, current, object$fixed,
        mesh = current$integration_mesh, information = FALSE)
      list(value = evaluate_model(theta, object$prepared, current, object$fixed),
        mesh = current$integration_mesh)
    }, error = identity)
    failed = inherits(result, 'error')
    if (name == 'person_effect' && native_outer && !failed &&
        identical(quadrature_outer_orders(result$value$diagnostics$native), outer_orders)) return(NULL)
    results[[name]] <<- if (!failed) result$value$diagnostics else NULL
    if (failed) failures[[name]] <<- list(class = class(result), message = conditionMessage(result),
      diagnostics = result[['diagnostics']], probes = result[['probes']])
    # Tightened targets can leave every cell of the fitted mesh in place; the
    # comparison is then a repeat evaluation and says so.
    changed = if (failed || is.null(fitted_mesh)) NA else
      !identical(result$mesh[c('outer', 'inner')], fitted_mesh[c('outer', 'inner')])
    comparison_row(name, current$nodes, control$nodes,
      if (failed) rep(NA_real_, 4L) else refinement_change(result$value, baseline, free, jacobian),
      failed, if (failed) conditionMessage(result) else if (isFALSE(changed))
        'refined targets left the fitted mesh unchanged' else '', tolerance,
      observation_nodes = if (inner) current$observation_nodes else NA_integer_,
      observation_quadrature = if (inner) current$observation_quadrature else NA_character_,
      mesh_changed = changed)
  })
  comparison = do.call(rbind, comparison)
  if (is.null(comparison)) stop('At least one effective quadrature order must increase', call. = FALSE)
  if ('both' %in% comparison$rule && !'person_effect' %in% comparison$rule) {
    comparison = comparison[comparison$rule != 'both', , drop = FALSE]
    if (is.null(nodes)) refined$nodes = control$nodes
  }
  refined$integration_mesh = refined$.rules = NULL
  # Refits go through update(): warm starts, subset rows and the call record.
  fitted = if (refit) do.call(stats::update, list(object, control = refined)) else NULL
  quadrature_check(object, 'higher_order', comparison, tolerance, fitted, integration = results,
    failures = failures,
    qualification = 'Same-method positive-mesh/order refinement at retained parameters; no accuracy certificate or silent refit.')
}

# GL with moving cuts: the fitted per-segment order becomes the reference.
check_split_quadrature = function(object, nodes, observation_nodes, tolerance, refit) {
  if (!is.null(observation_nodes))
    stop("GL with moving cuts has no PLN observation rule to refine", call. = FALSE)
  control = object$control
  refined = control
  refined$nodes = nodes %||% (control$nodes + max(5L, ceiling(control$nodes/2)))
  refined$quadrature_options$reference_nodes = control$nodes
  refined = validate_panel_control(refined)
  if (refined$nodes <= control$nodes)
    stop("The GL per-segment order must increase", call. = FALSE)
  theta = coef(object)
  baseline = evaluate_model(theta, object$prepared, control, object$fixed)
  if (length(baseline$diagnostics$dimension) && all(baseline$diagnostics$dimension == 0L))
    stop("This model has no numerical integral to refine", call. = FALSE)
  result = tryCatch(evaluate_model(theta, object$prepared, refined, object$fixed), error = identity)
  failed = inherits(result, "error")
  free = !object$parameter_info$fixed
  jacobian = parameter_jacobian(theta, object$parameter_info)[free]
  comparison = comparison_row("gl", refined$nodes, control$nodes,
    if (failed) rep(NA_real_, 4L) else refinement_change(result, baseline, free, jacobian),
    failed, if (failed) conditionMessage(result) else "", tolerance)
  # Refits go through update(): warm starts, subset rows and the call record.
  fitted = if (refit) do.call(stats::update, list(object, control = refined)) else NULL
  quadrature_check(object, "higher_order", comparison, tolerance, fitted,
    integration = if (!failed) result$diagnostics else NULL,
    qualification = "Finite GL-rule refinement; no integration-method substitution or global accuracy guarantee.")
}

# An independent GL rule at the fitted parameters of a fit whose person-effect
# rule is Gauss-Hermite, routed as panel_control(quadrature = "gl") routes the
# model (moving cuts or a fresh staged mesh); an active PLN observation integral
# also uses the GL rule. The fit's own rule is the reference row's baseline.
check_gl_reference = function(object, nodes, observation_nodes, tolerance, refit, options) {
  control = object$control; theta = coef(object); prepared = object$prepared
  if (uses_gl_quadrature(control))
    stop('compare = "gl": this fit already uses GL; use compare = "higher_order"', call. = FALSE)
  inner = prepared$family == "pln" && !is_fixed_zero("sd_error", theta, object$fixed) &&
    any(prepared$d == 1)
  if (model_dimension(theta, prepared, control, object$fixed) == 0L && !inner)
    stop("This model has no numerical integral to check", call. = FALSE)
  reference = unclass(control)
  reference[c("integration_mesh", ".rules")] = list(NULL)
  reference$quadrature = "gl"; reference$nodes = nodes %||% 20L
  reference$quadrature_options = options
  if (inner) {
    reference$observation_quadrature = "gl"; reference$observation_nodes = observation_nodes %||% 20L
    reference$observation_options = list()
  }
  # Invalid orders or quadrature_options stop; only integration may fail below.
  reference = validate_panel_control(structure(reference, class = "panelselect_control"))
  baseline = evaluate_model(theta, prepared, control, object$fixed)
  free = !object$parameter_info$fixed
  jacobian = parameter_jacobian(theta, object$parameter_info)[free]
  result = tryCatch({
    if (uses_staged_integration(prepared, reference, object$fixed)) {
      reference$integration_mesh = prepare_staged_integration(theta, prepared, reference,
        object$fixed, epoch = 1L)
      validate_staged_integration(theta, prepared, reference, object$fixed,
        mesh = reference$integration_mesh, information = FALSE)
    }
    evaluate_model(theta, prepared, reference, object$fixed)
  }, error = identity)
  failed = inherits(result, "error")
  comparison = comparison_row("gl", reference$nodes, control$nodes,
    if (failed) rep(NA_real_, 4L) else refinement_change(result, baseline, free, jacobian),
    failed, if (failed) conditionMessage(result) else "", tolerance,
    observation_nodes = if (inner) reference$observation_nodes else NA_integer_,
    observation_quadrature = if (inner) "gl" else NA_character_)
  reference$integration_mesh = reference$.rules = NULL
  fitted = if (refit && !failed) do.call(stats::update, list(object, control = reference)) else NULL
  quadrature_check(object, "gl", comparison, tolerance, fitted,
    integration = if (!failed) result$diagnostics else NULL,
    failures = if (failed) list(gl = list(class = class(result), message = conditionMessage(result),
      diagnostics = result[["diagnostics"]])) else list(),
    qualification = "Independent GL rule at the fitted parameters; no accuracy certificate or silent refit.")
}

# Changes from the fitted rule: total and largest panel log-likelihood change,
# largest total and largest panel optimizer-score change (free coordinates).
# 3. Shared comparison records and printing

refinement_change = function(result, baseline, free, jacobian) {
  delta = result$loglik - baseline$loglik
  score = sweep(result$score[, free, drop = FALSE] - baseline$score[, free, drop = FALSE],
    2L, jacobian, `*`)
  c(sum(delta), max(c(0, abs(delta))), max(c(0, abs(colSums(score)))), max(c(0, abs(score))))
}

# One comparison row with the columns shared by every route; NA marks a column
# that does not apply to the route.
comparison_row = function(rule, nodes, reference_nodes, change, failed, message, tolerance,
                          observation_nodes = NA_integer_, observation_quadrature = NA_character_,
                          mesh_changed = NA) {
  data.frame(rule = rule, nodes = paste(nodes, collapse = " x "),
    reference_nodes = paste(reference_nodes, collapse = " x "),
    observation_nodes = as.integer(observation_nodes),
    observation_quadrature = as.character(observation_quadrature),
    mesh_changed = as.logical(mesh_changed),
    loglik_change = change[1L], max_panel_change = change[2L], max_score_change = change[3L],
    max_panel_score_change = change[4L], integration_failed = failed,
    within_tolerance = if (failed) NA else all(abs(change[1:2]) <= tolerance[["loglik"]]) &&
      all(change[3:4] <= tolerance[["score"]]),
    message = message, row.names = NULL, stringsAsFactors = FALSE)
}

quadrature_check = function(object, compare, comparison, tolerance, refit, integration = NULL,
                            failures = list(), qualification) {
  structure(list(compare = compare, fitted = integration_description(object),
    comparison = comparison, tolerance = tolerance,
    optimizer_converged = object$optimizer$converged,
    covariance_valid = all(is.finite(stats::vcov(object))), integration = integration,
    failures = failures, refit = refit, qualification = qualification),
    class = "panelselect_quadrature_check")
}

#' @rdname check_quadrature
#' @param x A quadrature check to print.
#' @param digits Significant digits for the printed comparison; default `3`. The stored discrepancies retain full precision.
#' @export
print.panelselect_quadrature_check = function(x, digits = 3L, ...) {
  check = current_check(x)
  cat("Quadrature check (", check$compare, ") of: ", check$fitted, "\n", sep = "")
  shown = check$comparison[, c("rule", "nodes", "observation_nodes", "loglik_change",
    "max_panel_change", "max_score_change", "max_panel_score_change", "within_tolerance"),
    drop = FALSE]
  if (all(is.na(shown$observation_nodes))) shown$observation_nodes = NULL
  print(shown, digits = digits, row.names = FALSE)
  failed = check$comparison$integration_failed
  if (any(failed)) cat("Integration failed for:", paste(check$comparison$rule[failed], collapse = ", "),
    "(see $comparison$message)\n")
  cat("Tolerances: loglik ", format(check$tolerance[["loglik"]]), ", score ",
      format(check$tolerance[["score"]]),
      "; optimizer ", if (isTRUE(check$optimizer_converged)) "converged" else "did not report convergence",
      "; covariance ", if (isTRUE(check$covariance_valid)) "valid" else "unavailable",
      if (!is.null(check$refit)) "; refit stored in $refit", "\n", sep = "")
  invisible(x)
}

# Earlier saved checks name their comparison `method` (refine or gk), the
# rows outer, inner, joint and the moving-cut GL row split, and the observation
# columns inner_nodes and inner_quadrature, with the earlier rule values.
current_check = function(x) {
  if (is.null(x[["compare", exact = TRUE]]) && !is.null(x[["method", exact = TRUE]])) {
    x$compare = if (identical(x$method, "refine")) "higher_order" else x$method
    x$method = NULL
  }
  comparison = x$comparison
  rows = c(split = "gl", outer = "person_effect", inner = "observation", joint = "both")
  old_rows = comparison$rule %in% names(rows)
  comparison$rule[old_rows] = rows[comparison$rule[old_rows]]
  # Staged checks key their diagnostics and failures by row.
  for (part in intersect(c("integration", "failures"), names(x))) if (is.list(x[[part]])) {
    keyed = names(x[[part]]) %in% names(rows)
    names(x[[part]])[keyed] = rows[names(x[[part]])[keyed]]
  }
  earlier = c(inner_nodes = "observation_nodes", inner_quadrature = "observation_quadrature")
  hit = names(comparison) %in% names(earlier)
  names(comparison)[hit] = earlier[names(comparison)[hit]]
  rules = comparison$observation_quadrature
  if (is.character(rules)) {
    old = !is.na(rules) & rules %in% names(old_rule_values)
    comparison$observation_quadrature[old] = old_rule_values[rules[old]]
  }
  x$comparison = comparison
  x
}

# Each diagnostic leaf records the dimension and orders of its actual rule.
# Keep component paths: factorized outcomes can use the second supplied order.
quadrature_outer_orders = function(diagnostics) {
  if (!is.null(diagnostics$dimension))
    return(as.integer(rep_len(diagnostics$nodes, diagnostics$dimension)))
  unlist(lapply(Filter(is.list, diagnostics), quadrature_outer_orders), use.names = TRUE)
}

# 4. Stored stationarity diagnostics

#' Inspect score and curvature diagnostics without refitting
#'
#' Calculate stationarity diagnostics from the fitted object's stored scores and
#' covariances. The primary screen is `Q_OPG`; when valid observed information is
#' stored, `Q_Hessian` and local Newton diagnostics are also available. No
#' likelihood, score or Hessian is evaluated by this function.
#'
#' @param object A fitted `panelselect` model. An OPG-only fit normally has no
#'   observed information: assign `fit = add_information(fit)` if you explicitly
#'   want to calculate and retain it. Merely inspecting a fit does not add it.
#'
#' @section Diagnostic quantities:
#' Let \eqn{g} be the total free-parameter score in optimizer coordinates,
#' \eqn{B} the sum of independent-unit score outer products, and
#' \eqn{A=-\nabla^2\ell} the observed information for the fitted numerical objective.
#'
#' | Quantity | Definition and interpretation |
#' |---|---|
#' | `Q_OPG` | \eqn{g^\mathsf{T}B^{-1}g}. Screens score size relative to independent-unit score variation; it is not a Hessian-based distance. |
#' | `Q_Hessian` | \eqn{g^\mathsf{T}A^{-1}g}. Local curvature diagnostic, reported only with valid positive-definite information and inverse. |
#' | `newton_step` | \eqn{A^{-1}g}. Local quadratic estimate of the optimizer-coordinate displacement toward a stationary point. |
#' | `standardized_step` | Newton components divided by the square roots of the corresponding observed-inverse diagonals. |
#' | `max_abs_score`, `max_standardized_step` | Largest absolute score or standardized Newton component. |
#' | `information_eigenvalues`, `information_condition` | Curvature eigenvalues and, for positive-definite information, their maximum/minimum ratio. |
#'
#' Optimizer coordinates are coefficients, log-SDs and Fisher-transformed
#' correlations. The Hessian and Newton interpretation is local to these
#' coordinates; it is not generally invariant to nonlinear reparameterization
#' away from a stationary point. Neither Q is a p-value or a global-optimality test.
#'
#' @section Warning policy and unavailable results:
#' - The default warning cutoff is **`Q_OPG > 0.1`**, strictly. Change
#'   `panel_control(stationarity_warning = ...)` to use another Q cutoff;
#'   `NULL` disables the warning while preserving the statistic. This is a
#'   diagnostic policy, not an optimizer stopping criterion.
#' - There is no additional automatic cutoff for `Q_Hessian`.
#' - Missing, nonfinite, singular or indefinite observed information/inverse
#'   leaves `Q_Hessian` and the Newton quantities unavailable with a reason.
#'   OPG and sandwich covariance are never substituted for observed information.
#' - An all-fixed model has status `"all_fixed"` and no free-parameter
#'   stationarity quantities. An unavailable OPG screen is not an automatic pass.
#' - Native saved fits use current Q semantics and their stored scores/covariances.
#'   Previously stored Gaussian Hessians remain historical calculations; an
#'   explicit update/refit obtains corrected outcome-scale information.
#'
#' A small Q, or a backend convergence code, does not certify integration
#' accuracy, identification, calibrated intervals or a global maximum. Near
#' parameter boundaries the screen can overstate remaining displacement. Use
#' [check_quadrature()] for integration sensitivity, and explicit tighter
#' optimizer settings or different starts when investigating optimization.
#'
#' @return A list with top-level numeric `Q_OPG` and `Q_Hessian` (`NA` when
#'   unavailable), backend convergence/code/status, named optimizer `score`,
#'   the curvature/Newton quantities above, `information_error`, and retained
#'   rejected-trial summaries. `opg` contains `Q_OPG`, `status`, `tolerance`,
#'   `warning`, `reason`, and `basis`. Top-level `status` and `reason` concern
#'   the observed-information diagnostics; they can be unavailable while
#'   `opg$status` is available. The supplied fit is unchanged.
#' @family diagnostic functions
#' @seealso [panel_control()], [add_information()], [check_quadrature()]
#' @examples
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#' fit = panel_select(y ~ x, observed ~ x + z, dat, "person", "linear")
#' check_stationarity(fit)$opg
#' check_stationarity(add_information(fit))$max_standardized_step
#' @export
check_stationarity = function(object) {
  if (!inherits(object, "panelselect")) stop("object must be a panelselect fit")
  object = current_fit(object)
  optimizer = object[["optimizer"]]
  out = list(optimizer_converged = optimizer[["converged"]],
    optimizer_code = optimizer[["code"]], score = optimizer[["gradient"]],
    Q_OPG = NA_real_, Q_Hessian = NA_real_,
    max_abs_score = NA_real_, information_eigenvalues = NULL,
    information_condition = NA_real_, newton_step = NULL,
    standardized_step = NULL, max_standardized_step = NA_real_,
    information_error = optimizer[["information_error"]],
    opg = NULL, rejections = optimizer[["rejections"]],
    optimizer_status = optimizer[["status"]], status = "unavailable", reason = NULL)
  unavailable = function(reason) { out$reason = reason; out }
  info = object[["parameter_info"]]
  if (!is.data.frame(info) || !all(c("name", "fixed") %in% names(info)))
    return(unavailable("Free-parameter metadata is unavailable"))
  parameters = info$name[!info$fixed]
  p = length(parameters)
  if (!p) {
    out$opg = opg_stationarity(numeric(), matrix(numeric(), 0L, 0L))
    out$status = "all_fixed"; out$reason = "There are no free optimizer parameters"
    return(out)
  }
  score = out$score
  if (!is.numeric(score) || length(score) != p || any(!is.finite(score)) ||
      anyDuplicated(names(score)) || !setequal(names(score), parameters))
    return(unavailable("A finite, named free-parameter score is unavailable"))
  out$score = score[parameters]
  out$max_abs_score = max(abs(out$score))
  tolerance = if ("stationarity_warning" %in% names(object$control))
    object$control$stationarity_warning else 0.1
  out$opg = opg_stationarity(out$score, optimizer[["covariance"]][["opg"]], tolerance)
  out$Q_OPG = out$opg$Q_OPG
  if (identical(optimizer[["information_computed"]], FALSE))
    return(unavailable("Observed information was not computed"))
  hessian = hessian_stationarity(out$score, optimizer[["information"]],
    optimizer[["covariance"]][["observed"]])
  out[names(hessian)] = hessian
  out
}

# Stored curvature diagnostics share one validation path with fitted records.
# Both the score and information use the same free optimizer coordinates.
hessian_stationarity = function(score, information, covariance) {
  out = list(Q_Hessian = NA_real_, information_eigenvalues = NULL,
    information_condition = NA_real_, newton_step = NULL,
    standardized_step = NULL, max_standardized_step = NA_real_,
    status = "unavailable", reason = NULL)
  unavailable = function(reason) { out$reason = reason; out }
  parameters = names(score)
  p = length(score)
  if (!p) {
    out$status = "all_fixed"
    return(unavailable("There are no free optimizer parameters"))
  }
  if (!is.numeric(score) || is.complex(score) || any(!is.finite(score)) ||
      is.null(parameters) || anyNA(parameters) || anyDuplicated(parameters))
    return(unavailable("A finite, named free-parameter score is unavailable"))
  aligned_matrix = function(value) is.matrix(value) && is.numeric(value) &&
    !is.complex(value) && identical(dim(value), c(p, p)) && all(is.finite(value)) &&
    !anyDuplicated(rownames(value)) && !anyDuplicated(colnames(value)) &&
    setequal(rownames(value), parameters) && setequal(colnames(value), parameters)
  if (!aligned_matrix(information))
    return(unavailable("Finite, named observed information is unavailable"))
  information = information[parameters, parameters, drop = FALSE]
  if (!isSymmetric(information))
    return(unavailable("Stored observed information is not symmetric"))
  eigenvalues = eigen(information, symmetric = TRUE, only.values = TRUE)$values
  out$information_eigenvalues = eigenvalues
  if (min(eigenvalues) <= 0)
    return(unavailable("Observed information is not positive definite"))
  out$information_condition = max(eigenvalues)/min(eigenvalues)
  if (!aligned_matrix(covariance))
    return(unavailable("A finite, named observed-information inverse is unavailable"))
  covariance = covariance[parameters, parameters, drop = FALSE]
  if (!isSymmetric(covariance) || any(diag(covariance) <= 0))
    return(unavailable("Stored observed-information inverse is invalid"))
  root = tryCatch(chol(covariance), error = function(error) NULL)
  if (is.null(root))
    return(unavailable("Stored observed-information inverse is not positive definite"))
  quadratic = sum((root %*% score)^2)
  if (!is.finite(quadratic))
    return(unavailable("The observed-information quadratic is not finite"))
  step = drop(covariance %*% score)
  standardized = step/sqrt(diag(covariance))
  if (any(!is.finite(standardized)))
    return(unavailable("The local Newton displacement is not finite"))
  out$newton_step = step
  out$standardized_step = standardized
  out$max_standardized_step = max(abs(standardized))
  out$Q_Hessian = quadratic
  out$status = "available"
  out
}

# This metric screens the total score using independent-panel score variation.
# It is not an observed-Hessian Newton step or an optimizer stopping condition.
opg_stationarity = function(gradient, covariance, tolerance = 0.1) {
  out = list(status = "unavailable", Q_OPG = NA_real_,
    tolerance = tolerance, warning = FALSE, reason = NULL, basis = "panel_opg")
  unavailable = function(reason) { out$reason = reason; out }
  p = length(gradient)
  if (!p) {out$status = "all_fixed"; return(out)}
  if (!is.numeric(gradient) || is.complex(gradient) || any(!is.finite(gradient)) || is.null(names(gradient)) ||
      anyDuplicated(names(gradient))) return(unavailable("A finite named optimizer score is unavailable"))
  if (!is.matrix(covariance) || !is.numeric(covariance) || is.complex(covariance) ||
      !identical(dim(covariance), c(p, p)) || any(!is.finite(covariance)) ||
      anyDuplicated(rownames(covariance)) || anyDuplicated(colnames(covariance)) ||
      !setequal(rownames(covariance), names(gradient)) ||
      !setequal(colnames(covariance), names(gradient)))
    return(unavailable("A finite named OPG inverse is unavailable"))
  covariance = covariance[names(gradient), names(gradient), drop = FALSE]
  if (!isSymmetric(covariance)) return(unavailable("The OPG inverse is not symmetric"))
  root = tryCatch(chol(covariance), error = function(error) NULL)
  if (is.null(root)) return(unavailable("The OPG inverse is not positive definite"))
  quadratic = sum((root %*% gradient)^2)
  if (!is.finite(quadratic)) return(unavailable("The OPG quadratic is not finite"))
  out$Q_OPG = quadratic
  out$status = "available"
  out$warning = !is.null(tolerance) && out$Q_OPG > tolerance
  out
}

# Free correlations near +-1 and random-effect SDs near zero lie where the atanh
# and log coordinates drift along flat ridges, so the OPG screen can overstate
# the remaining distance to an optimum. A Linear outcome
# SD is judged relative to the residual SD; other latent scales are unit free.
# Returns NULL when no free parameter is near its boundary.
boundary_note = function(theta, info, model) {
  free = !info$fixed
  reference = rep(1, length(theta))
  if ("sd_error" %in% model$positive_scales && "sd_re_outcome" %in% names(theta))
    reference[names(theta) == "sd_re_outcome"] = theta[["sd_error"]]
  near = free & ((info$type == "correlation" & abs(theta) > 0.99) |
    (info$type == "sd" & !names(theta) %in% model$positive_scales & theta < 1e-3*reference))
  if (!any(near)) return(NULL)
  paste0(paste0(names(theta)[near], " = ", signif(theta[near], 6), collapse = ", "),
    if (sum(near) == 1L) " is" else " are",
    " near the parameter boundary, where this screen can overstate the remaining distance")
}
