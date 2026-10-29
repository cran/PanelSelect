# Summaries and printed views of a fit: coefficient tables with the exclusion
# restriction, covariance, integration and optimizer descriptions.

# 1. Public summary and print methods

#' Summarize and print a fitted model
#'
#' Display the model specification, retained sample, estimates, uncertainty and
#' numerical diagnostics. `print(fit)` gives a compact coefficient view;
#' `summary(fit)` constructs the full coefficient table and its printed summary.
#'
#' @param object A fitted `panelselect` model to summarize.
#' @param x A fitted model for `print.panelselect()`, or the result of
#'   `summary()` for `print.summary.panelselect()`.
#' @param vcov Covariance used by the summary: `"opg"`, `"observed"`,
#'   `"sandwich"`, or a named natural-scale covariance matrix. Defaults to the
#'   choice stored in the fit. Observed/sandwich requests calculate information
#'   if needed and reuse that calculation for `Q_Hessian`; the input fit is not
#'   modified. Assign [add_information()] first to retain it for later requests.
#' @param level Confidence level between zero and one for returned coefficient
#'   intervals; default `0.95`. Printed coefficient tables show estimates, SEs,
#'   z statistics and p-values; interval columns remain in the summary object.
#' @param digits Significant digits for printing; default is at least three,
#'   otherwise `getOption("digits") - 3`.
#' @param signif.stars Whether to print significance stars; defaults to
#'   `getOption("show.signif.stars")`.
#' @param ... `summary()` rejects unused arguments. Print methods ignore them;
#'   select covariance in `summary()`, not in `print()`.
#'
#' @section What to inspect:
#' | Component | Interpretation |
#' |---|---|
#' | Specification and sample | Outcome/selection formulas, panel versus cross-sectional structure, retained rows/people and omitted rows. |
#' | Exclusion restriction | Raw selection variables absent from the active outcome equation; not a statistical test of validity. |
#' | Coefficient table | Natural-scale estimates and SEs; correlations use Fisher-scale tests, SDs have no ordinary zero test, and fixed parameters have no inferential entries. |
#' | Structural-parameter legend | Meanings of the SD and correlation parameters for this model. |
#' | Covariance | Requested estimator and its availability; it is never silently replaced. See [panelselect-inference]. |
#' | Integration | The fitted rules/orders and available integration diagnostics. See [check_quadrature()]. |
#' | Optimizer and Q | Backend stopping status, `Q_OPG`, and `Q_Hessian` when valid information exists. These assess different aspects of the fitted numerical objective. |
#'
#' `Q_OPG > stationarity_warning` is flagged as substantial (`0.1` by default).
#' There is no automatic warning cutoff for `Q_Hessian`. Neither a small Q nor
#' optimizer convergence proves integration accuracy, identification or a global
#' maximum. Saved native fits are displayed under the current Q convention using
#' stored scores/covariances; inspecting them does not refit or overwrite them.
#'
#' @return `summary()` returns a `summary.panelselect` list. `coefficients`
#'   contains `estimate`, `se`, `z`, `p`, `lci`, and `uci`; `stationarity` contains
#'   the OPG status/threshold/flag and the two named Q values. Other components
#'   describe the model, likelihood criteria, covariance and numerical settings.
#'   `exclusion_restriction` is a character vector (`character(0)` when none,
#'   `NULL` for standalone models). Print methods return their input invisibly.
#' @examples
#' dat = data.frame(person = rep(1:20, each = 3), y = rep(c(0, 1, 3, 1, 2, 4), 10))
#' fit = panel_select(y ~ 1, NULL, dat, "person", "poisson",
#'                    fixed = c(sd_re_outcome = 0))
#' print(fit)
#' answer = summary(fit, vcov = "observed", level = 0.9)
#' print(answer, digits = 4, signif.stars = FALSE)
#' answer$stationarity
#' @rdname summary.panelselect
#' @seealso [panelselect-methods], [check_stationarity()], [panelselect-inference]
#' @export
summary.panelselect = function(object, vcov = object$control$vcov, level = 0.95, ...) {
  object = current_fit(object)
  reject_unused(list(...), "summary",
    "the covariance is selected with vcov = and the interval level with level =")
  type = if (is.matrix(vcov)) "supplied" else match.arg(vcov, c("observed", "opg", "sandwich"))
  # Keep requested information locally so inference and Q_Hessian share it.
  # The caller's fit is unchanged; OPG/supplied covariance needs no new Hessian.
  if (type %in% c("observed", "sandwich") &&
      identical(object$optimizer[["information_computed"]], FALSE))
    object = add_information(object)
  covariance = covariance_request(object, vcov)
  diagnostic = check_stationarity(object)
  stationarity = diagnostic$opg
  stationarity$Q_Hessian = diagnostic$Q_Hessian
  free = !object$parameter_info$fixed
  available = all(is.finite(covariance$natural[free, free, drop = FALSE]))
  structure(c(list(call = object$call, model = object$prepared$model,
                   model_id = object[["model_id"]] %||% object$prepared[["model_id"]],
                   coefficients = coefficient_table(object, covariance, level),
                   parameter_legend = parameter_legend(object),
                   fixed = names(object$fixed), covariance = type,
                   covariance_status = if (available) "available" else covariance_problem(object, type),
                   loglik = object$loglik, df = sum(free),
                   AIC = AIC.panelselect(object), BIC = BIC.panelselect(object),
                   converged = object$optimizer$converged,
                   optimizer_status = object$optimizer[["status"]],
                   stationarity = stationarity,
                   integration = object$integration, time = object$time),
                fit_description(object, diagnostic)),
            class = "summary.panelselect")
}

#' @rdname summary.panelselect
#' @export
print.panelselect = function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  x = current_fit(x)
  description = fit_description(x)
  print_fit_header(description)
  free = !x$parameter_info$fixed
  cat("\nCoefficients:\n")
  if (any(free)) print(x$coefficients[free], digits = digits) else cat("No free parameters.\n")
  print_fit_fixed(description)
  cat("Log likelihood: ", format(x$loglik, digits = digits), "\n", sep = "")
  print_fit_status(description)
  invisible(x)
}

#' @rdname summary.panelselect
#' @export
print.summary.panelselect = function(x, digits = max(3L, getOption("digits") - 3L),
                                      signif.stars = getOption("show.signif.stars"), ...) {
  # Saved summaries contain parameter labels but no fit to normalize. Change
  # those labels locally; keep formulas, regressor names and the ID untouched.
  rownames(x$coefficients) = current_names(rownames(x$coefficients))
  x$fixed = current_names(x$fixed)
  if (length(x$parameter_legend)) names(x$parameter_legend) = current_names(names(x$parameter_legend))
  for (old in names(renamed_parameters)) if (length(x$fixed_values))
    x$fixed_values = gsub(paste0("(?<![[:alnum:]_.])", old, "(?![[:alnum:]_.])"),
      renamed_parameters[[old]], x$fixed_values, perl = TRUE)
  # A summary saved before these fields existed prints its stored essentials.
  if (is.null(x$label)) x$label = paste("PanelSelect", x[["model_id"]] %||% x$model, "model")
  # Old summaries do not retain the scores needed by check_stationarity().
  # Their saved quadratic suffices to display the new Q convention.
  screen = x$stationarity
  if (any(c("score_norm", "quadratic") %in% names(screen))) {
    screen$Q_OPG = screen$quadratic %||% screen$score_norm^2
    if (!is.finite(screen$Q_OPG)) screen$Q_OPG = NA_real_
    screen[c("score_norm", "quadratic")] = NULL
    screen$Q_Hessian = NA_real_
    screen$warning = is.finite(screen$Q_OPG) && !is.null(screen$tolerance) &&
      screen$Q_OPG > screen$tolerance
    x$stationarity = screen
    x$optimizer_description = sub("; OPG-scaled score .*", "", x$optimizer_description)
    if (is.finite(screen$Q_OPG))
      x$optimizer_description = paste0(x$optimizer_description, "; Q_OPG = ",
        format(signif(screen$Q_OPG, 4)), if (screen$warning) " (substantial)")
  }
  print_fit_header(x)
  table = x$coefficients
  free = !rownames(table) %in% x$fixed
  cat("\n")
  if (any(free)) {
    stats::printCoefmat(table[free, 1:4, drop = FALSE], digits = digits,
                        signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "")
    rows = rownames(table)[free]
    notes = c(if (any(startsWith(rows, "cor_"))) "Correlation z and p values use the Fisher (atanh) scale",
              if (any(startsWith(rows, "sd_"))) "SD rows have no zero test")
    if (length(notes)) cat(paste(notes, collapse = "; "), ".\n", sep = "")
  } else cat("No free parameters.\n")
  legend = x[["parameter_legend"]]
  if (length(legend)) cat("\nStructural parameters:\n",
    paste0("  ", format(names(legend)), "  ", legend, "\n"), sep = "")
  print_fit_fixed(x)
  cross = identical(x[["structure"]], "cross_sectional")
  cat("\nLog likelihood: ", format(x$loglik, digits = digits),
      if (!is.null(x$df)) paste0(" (df ", x$df, ")"),
      "; AIC ", format(x$AIC, digits = digits),
      "; BIC (", if (cross) "observations" else "panels", ") ",
      format(x$BIC, digits = digits), "\n", sep = "")
  cat("Covariance:  ", covariance_label(x$covariance, cross),
      if (!is.null(x$covariance_status)) paste0(", ", x$covariance_status), "\n", sep = "")
  print_fit_status(x)
  invisible(x)
}

# 2. Shared fit descriptions

# Descriptive fields shared by print() and summary(); no covariance is computed.
fit_description = function(object, diagnostic = check_stationarity(object)) {
  prepared = object$prepared
  structure = model_structure(object)
  list(label = paste("PanelSelect fit:", model_label(object[["model_id"]] %||% prepared[["model_id"]])),
       structure = structure,
       formulas = list(outcome = formula.panelselect(object),
                       selection = if (prepared$has_selection) formula.panelselect(object, "selection")),
       id = if (structure == "panel") prepared$id,
       n_rows = prepared$n_rows, n_panels = prepared$n_panels,
       omitted = length(prepared$na.action),
       exclusion_restriction = exclusion_restriction(object),
       fixed_values = fixed_description(object),
       integration_description = integration_description(object),
       optimizer_description = optimizer_description(object, diagnostic))
}

print_fit_header = function(x) {
  cat(x$label, "\n", sep = "")
  if (!is.null(x$formulas$outcome)) cat("Outcome:   ", deparse1(x$formulas$outcome), "\n", sep = "")
  if (!is.null(x$formulas$selection)) cat("Selection: ", deparse1(x$formulas$selection), "\n", sep = "")
  sample = if (identical(x[["structure"]], "cross_sectional"))
    paste(x$n_rows, "independent observations") else
    paste0(x$n_rows, " rows in ", x$n_panels, " panels", if (!is.null(x$id)) paste0(" (id \"", x$id, "\")"))
  omitted = x$omitted %||% 0L
  cat("Sample:    ", sample, if (omitted) paste0("; ", omitted, " row", if (omitted > 1L) "s",
    " omitted (missing data)"), "\n", sep = "")
  restriction = x[["exclusion_restriction"]]
  if (!is.null(restriction)) cat("Exclusion restriction: ", if (length(restriction))
    paste(restriction, collapse = ", ") else
      "none (every selection variable also enters the outcome equation)", "\n", sep = "")
}

print_fit_fixed = function(x)
  if (length(x$fixed_values)) cat("Fixed: ", paste(x$fixed_values, collapse = "; "), "\n", sep = "")

print_fit_status = function(x) {
  if (!is.null(x$integration_description))
    cat("Integration: ", x$integration_description, "\n", sep = "")
  if (!is.null(x$optimizer_description))
    cat("Optimizer:   ", x$optimizer_description, "\n", sep = "")
}

# 3. Model and inference labels

# What each structural parameter of the model means, in parameter order.
parameter_legend = function(object) {
  family = object$prepared$family
  outcome = switch(family, gaussian = "outcome error", probit = "outcome probit error",
    "row-level normal error in the log outcome rate")
  text = c(sd_re_selection = "SD of the person effect in the selection equation",
    sd_re_outcome = "SD of the person effect in the outcome equation",
    sd_error = paste("SD of the", outcome),
    cor_re = "correlation of the selection and outcome person effects",
    cor_error = paste("correlation of the selection error and the", outcome, "in the same row"))
  text[intersect(object$parameter_info$name, names(text))]
}

# Readable model names; an ASCII arrow separates the selection and outcome equations.
model_label = function(model_id) {
  labels = c(probit_re = "Probit RE (standalone panel)",
    poisson_re = "Poisson RE (standalone panel)", pln_re = "PLN RE (standalone panel)",
    probit_re_probit_re = "Probit RE -> Probit RE (panel selection)",
    probit_re_linear_re = "Probit RE -> Linear RE (panel selection)",
    probit_re_poisson_re = "Probit RE -> Poisson RE (panel selection)",
    probit_re_pln_re = "Probit RE -> PLN RE (panel selection)",
    probit_probit = "Probit -> Probit (cross-sectional selection)",
    probit_linear = "Probit -> Linear (cross-sectional selection)",
    probit_pln = "Probit -> PLN (cross-sectional selection)")
  if (!is.null(model_id) && model_id %in% names(labels)) unname(labels[[model_id]]) else
    paste(model_id, "model")
}

# Selection variables absent from the outcome equation: character(0) when none,
# NULL without a selection equation. A term enters an equation when one of its
# columns has a coefficient not fixed at zero; offsets always enter. Variables
# are raw data variables, so poly(x, 2) and x share x.
exclusion_restriction = function(object) {
  prepared = object$prepared
  if (!prepared$has_selection) return(NULL)
  zero = names(object$fixed)[object$fixed == 0]
  variables = function(equation) {
    design = prepared$design[[equation]]
    # Fits saved before designs stored `assign` recover it from the model matrix.
    assign = design$assign %||% attr(model.matrix.panelselect(object, equation), "assign")
    free = !paste0(equation, ".", design$columns) %in% zero
    labels = attr(design$terms, "term.labels")[setdiff(unique(assign[free]), 0L)]
    offsets = as.list(attr(design$terms, "variables"))[-1L][attr(design$terms, "offset")]
    as.character(unique(c(unlist(lapply(labels, function(label) all.vars(str2lang(label)))),
                          unlist(lapply(offsets, all.vars)))))
  }
  setdiff(variables("selection"), variables("outcome"))
}

# Fixed values in parameter order; a zero correlation forced by a zero scale
# says so.
fixed_description = function(object) {
  fixed = object$fixed
  if (!length(fixed)) return(character())
  fixed = fixed[intersect(object$parameter_info$name, names(fixed))]
  prepared = object$prepared
  model = resolve_model(prepared$family, prepared$has_selection, model_structure(object),
    object[["model_id"]] %||% prepared[["model_id"]])
  zero = names(fixed)[fixed == 0]
  vapply(names(fixed), function(name) {
    scales = intersect(model$zero_correlations[[name]], zero)
    paste0(name, " = ", format(fixed[[name]]),
      if (length(scales)) paste0(" (required by ", paste0(scales, " = 0", collapse = ", "), ")"))
  }, character(1), USE.NAMES = FALSE)
}

covariance_label = function(type, cross = FALSE) {
  unit = if (cross) "observation" else "panel"
  switch(type %||% "",
    opg = paste0("OPG (inverse sum of ", unit, "-score outer products)"),
    observed = "inverse observed information (negative log-likelihood Hessian)",
    sandwich = paste0("sandwich (observed information with ", unit, "-score meat)"),
    supplied = "supplied natural-scale matrix",
    type %||% "")
}

# Why a covariance estimate is unavailable; it is never replaced by another.
covariance_problem = function(object, type) {
  error = object$optimizer[["information_error"]]
  reason = if (type %in% c("observed", "sandwich") && !is.null(error))
    paste("observed information failed:", error) else switch(type,
      opg = "the outer product of scores is singular",
      supplied = "the supplied matrix has non-finite free entries",
      "the observed information is singular, not positive definite or not computable")
  paste("unavailable:", reason)
}

# The integration actually used: the outer random-effect rule (none when no
# latent effect remains), the PLN observation rule and staged mesh use.
# 4. Numerical integration and optimizer status

integration_description = function(object) {
  control = object$control; integration = object$integration
  theta = object$coefficients; prepared = object$prepared
  q = model_dimension(theta, prepared, control, object$fixed)
  outer = if (q == 0L) NULL else if (identical(integration$method, "split_legendre"))
    paste("Gauss-Legendre (gl) with moving cuts, order", integration$canonical_order,
      "per segment; refinement order", integration$reference_order) else
    if (identical(control$quadrature, "gk")) removed_gk_label else
    if (identical(control$quadrature, "gl")) paste("staged Gauss-Legendre (gl), order", control$nodes) else {
      orders = rep_len(control$nodes, 2L)
      effects = if (q == 2L) c("selection", "outcome") else single_effect(theta, prepared)
      method = if (!is.null(control$fixed_geometry)) "frozen adaptive Gauss-Hermite (removed)" else
        if (control$quadrature == "agh") "adaptive Gauss-Hermite (agh)" else "fixed Gauss-Hermite (gh)"
      paste0(method, ", ",
        paste(orders[match(effects, c("selection", "outcome"))], collapse = " x "), " nodes (",
        if (q == 2L) "selection x outcome" else paste(effects, "effect"), ")")
    }
  observation = if (identical(prepared$family, "pln") && !is_fixed_zero("sd_error", theta, object$fixed))
    switch(control$observation_quadrature %||% "gh",
      gh = paste0("fixed Gauss-Hermite (gh), ", control$observation_nodes, " nodes"),
      agh = paste0("adaptive Gauss-Hermite (agh), ", control$observation_nodes, " nodes"),
      gl = paste("staged Gauss-Legendre (gl), order", control$observation_nodes),
      gk = removed_gk_label)
  if (is.null(outer) && is.null(observation)) return("None needed (closed-form likelihood)")
  stages = integration$staging
  text = paste0(outer %||% "no random-effect integral",
    if (!is.null(observation)) paste("; PLN observation rule:", observation),
    if (!is.null(stages)) paste0("; ", stages$stages_used, " mesh stage",
      if (stages$stages_used != 1L) "s", " used of at most ", stages$maximum_stages))
  paste0(toupper(substring(text, 1L, 1L)), substring(text, 2L))
}

# Backend status, work and the separate OPG stationarity screen.
optimizer_description = function(object, diagnostic = check_stationarity(object)) {
  optimizer = object$optimizer; control = object$control
  status = optimizer[["status"]] %||% if (isTRUE(optimizer$converged)) "converged" else "not_converged"
  # optim reports function and gradient evaluations; maxLik reports iterations.
  counts = optimizer$counts
  work = if (identical(control$optimizer, "maxLik") && length(counts)) paste(counts[[1L]], "iterations") else
    if (length(counts) >= 2L) paste(counts[[1L]], "function and", counts[[2L]], "gradient evaluations") else
      "evaluation counts unavailable"
  stages = object$integration$staging$stages_used
  if (!is.null(stages) && stages > 1L && !is.null(optimizer$evaluations))
    work = paste0(work, " in the final mesh stage; ", optimizer$evaluations,
      " likelihood evaluations over ", stages, " stages")
  backend = paste(control$optimizer, control$method)
  text = switch(status,
    converged = paste0("converged (", backend, "; ", work, ")"),
    not_converged = paste0("did not report convergence (", backend, ", code ", optimizer$code, "; ", work, ")"),
    not_run = "not run (maxit = 0)",
    all_fixed = "not needed (all parameters fixed)",
    no_convergence_test = paste0("stopped (", backend, " has no convergence test; ", work, ")"),
    status)
  # Recompute flags from stored scores/covariances, including for saved fits.
  screen = diagnostic$opg
  if (identical(screen$status, "available")) {
    note = if (isTRUE(screen$warning)) boundary_note(object$coefficients, object$parameter_info,
      resolve_model(object$prepared$family, object$prepared$has_selection, model_structure(object),
        object[["model_id"]] %||% object$prepared[["model_id"]]))
    text = paste0(text, "; Q_OPG = ", format(signif(screen$Q_OPG, 4)),
      if (isTRUE(screen$warning)) paste0(" (substantial", if (!is.null(note)) paste0("; ", note), ")"))
  }
  if (identical(screen$status, "unavailable"))
    text = paste0(text, "; OPG stationarity screen unavailable: ", screen$reason)
  if (is.finite(diagnostic$Q_Hessian))
    text = paste0(text, "; Q_Hessian = ", format(signif(diagnostic$Q_Hessian, 4)))
  text
}
