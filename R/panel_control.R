# Numerical settings of a fit: panel_control(), its print method and option
# validators, the route predicates, the effective fit control, and the
# normalization and revalidation of stored controls.

# 1. Public controls and their printed description

#' Control numerical integration and optimization
#'
#' Collects numerical settings for a panel selection fit with [panel_select()],
#' also used by its cross-sectional selection and standalone models. Person effects
#' use fixed Gauss-Hermite (GH), 10 nodes per effect, by default; the separate
#' PLN observation integral uses fixed GH with 20 nodes. The default optimizer
#' is `optim` BFGS with 500 iterations. These are starting settings to assess,
#' not guarantees of accuracy. Missing-data handling, initialization and the
#' reported covariance belong to [panel_select()] (`na_action`, `init`, `vcov`).
#'
#' @param quadrature Person-effect integration rule: `"gh"` (fixed GH, default),
#'   `"agh"` (adaptive GH), or `"gl"` (Gauss-Legendre on segments). See
#'   Integration rules for how they differ. Cross-sectional models have no
#'   person-effect integral.
#' @param nodes Positive integer GH order per person effect, or one GL order
#'   per segment between 3 and 256. `NULL` chooses 10 for GH/AGH and 20 for GL.
#'   GH/AGH also accept two orders in selection/outcome order, optionally named
#'   `c(selection = 8, outcome = 12)`; GL requires one order.
#' @param quadrature_options Uniquely named list of GL person-effect options;
#'   see GL options for every name, default and meaning. Must be empty for
#'   GH/AGH. A supplied list overrides only its named GL option defaults.
#' @param observation_quadrature PLN observation-effect rule: `"gh"` (default),
#'   `"agh"` or `"gl"`, selected independently of `quadrature`. Other families
#'   do not use this rule; controls are still validated.
#' @param observation_nodes Positive integer PLN observation order. `NULL`
#'   chooses 20 for every rule; for GL it is the per-segment order, from 3 to 256.
#' @param observation_options Uniquely named list of GL observation options.
#'   It accepts the same names as `quadrature_options`, with a smaller default
#'   `batch_nodes`. Must be empty for a GH/AGH observation rule.
#' @param mode_maxit Positive integer iteration limit for adaptive-GH person-effect
#'   mode searches (50). It does not control the model optimizer or the separate
#'   PLN observation-effect mode search.
#' @param mode_tol Positive residual tolerance of adaptive-GH person-effect mode searches
#'   (1e-9). Non-PLN person-effect modes are polished to `min(mode_tol, 1e-12)`,
#'   accepting a roundoff floor within `mode_tol`; PLN modes use `mode_tol`.
#'   A relative residual above `10 * mode_tol` fails the evaluation.
#' @param exact_reductions Logical, default `TRUE`. Use exact likelihood and
#'   integration reductions implied by the model and fixed parameters. `FALSE`
#'   requests the general route for numerical comparisons; it does not change
#'   the restrictions or restore an unidentified parameter.
#' @param optimizer Backend: `"optim"` (default) or `"maxLik"`. The latter
#'   requires the optional maxLik package. See Optimizer methods.
#' @param method Backend method name, default `"BFGS"`; the available names
#'   differ by backend and are listed below.
#' @param maxit Nonnegative integer optimizer iteration limit (500), passed as
#'   `maxit` to `optim` and `iterlim` to maxLik. Zero evaluates the starting
#'   point without claiming optimizer convergence.
#' @param reltol Positive relative stopping tolerance (default
#'   `sqrt(.Machine$double.eps)`). Passed to maxLik and to `optim` methods other
#'   than `L-BFGS-B`; the automatic fixed-GH binary BFGS exception is below.
#' @param factr Positive `optim` `L-BFGS-B` function-reduction tolerance factor
#'   (1e7), relative to machine precision. Other supported methods do not use it.
#' @param parscale `NULL` for automatic BFGS scaling, `1` for an unscaled metric,
#'   or positive finite optimizer-coordinate scales. Supply one value or one
#'   per free parameter; a named vector must match the free parameter names.
#'   Native maxLik Newton/BHHH/BFGSR methods do not accept `parscale`.
#' @param trace Nonnegative integer backend tracing level; 0 is silent.
#'   Passed as `trace` to `optim` and `printLevel` to maxLik.
#' @param stationarity_warning Positive threshold for `Q_OPG = g' B^{-1} g`,
#'   where `g` is the total free-parameter score and `B` its summed
#'   independent-unit score outer products in optimizer coordinates. The
#'   default warns strictly when `Q_OPG > 0.1`. `NULL` disables this warning,
#'   while retaining the diagnostic. This is separate from backend stopping,
#'   integration checks and covariance validity; see [check_stationarity()].
#' @param ... Not used. Removed settings (`binary_integration`, `checked`,
#'   `fixed_geometry`, and GK), misplaced fit settings and earlier argument
#'   names stop with messages naming their replacements. Other names are errors.
#' @param x A `panelselect_control` list to print. Printing shows the selected
#'   integration rules, optimizer settings and Q warning threshold.
#'
#' @section Integration rules:
#' | Rule | Calculation | When to compare it |
#' |---|---|---|
#' | `"gh"` (default) | Fixed standard-normal nodes; order stays fixed during optimization. | Establish an initial fit, reproduce a fixed-GH specification, and assess order sensitivity. |
#' | `"agh"` | Centres and scales nodes at each panel's posterior mode and curvature; scores include node, curvature and volume derivatives. | Assess whether concentrating nodes around the posterior improves a difficult GH approximation. |
#' | `"gl"` | GL on segments, with a distinct order for empirical refinement checks and separate analytic tail envelopes. | Cross-check GH/AGH at a common parameter vector or fit with an explicitly selected GL rule. |
#'
#' Analytical scores differentiate the finite rule actually optimized. Correct
#' differentiation and backend convergence do not establish accuracy of the
#' continuous integral. Higher order is not uniformly more accurate; inspect
#' likelihood and score sensitivity with [check_quadrature()], especially for
#' large scales, near-boundary correlations and rare response histories.
#' Equal GH/GL orders are neither equal work nor matched accuracy, because GL
#' uses that many nodes in each segment.
#'
#' PLN has two integration levels. Choosing outer AGH or GL does not change the
#' fixed-GH observation default. Coarse fixed observation nodes can make the
#' person-effect integrand corrugated; `observation_quadrature = "agh"`
#' differentiates its own moving observation geometry and can help in some
#' cases, with additional computation. Select both rules explicitly for GL at
#' both levels. Refine observation and person-effect orders separately and then
#' jointly before interpreting changed estimates.
#'
#' Exact reductions skip integrals made unnecessary by fixed zero scales or
#' correlations. Cross-sectional selection without PLN is closed form. A
#' selected rule can therefore do no integration work in a reduced model.
#'
#' @section GL fitting and limits:
#' Standalone panel probit and selected binary/Linear panels with `cor_error`
#' fixed at zero use GL with analytically differentiated moving cuts. Other GL
#' routes keep a positive mesh fixed within each optimization and information
#' stage. A failed staged gate at any trial point ends that stage; the next
#' stage starts from its best valid point on a refined mesh. At most three
#' stages are allowed, including stages needed for requested covariance.
#' Unresolved preparation or exhausted limits fail explicitly and retain
#' diagnostic records. No rule, tolerance or covariance estimator is silently
#' substituted.
#'
#' Staged GL can be substantially more expensive than GH/AGH because it
#' prepares, checks and refines meshes. Cost depends on the model, active
#' dimensions, restrictions and accuracy target; the moving-cut routes and
#' staged PLN routes do not share a universal timing ratio. Physical moving-cut
#' GL can also be difficult when the person-effect correlation approaches -1
#' or 1. Refinement agreement remains empirical evidence, not a certified
#' integration-error bound.
#'
#' GL observed information must pass step/order checks and, when positive
#' definite, inverse/SE stability checks; otherwise it remains unavailable.
#' GL supports population prediction, serialization and explicit update/refit.
#' History-conditioned prediction and [random_effects()] are unavailable when
#' GL is selected for person effects or the PLN observation rule. There is no
#' prediction-time integration override; use an explicit GH/AGH refit if those
#' targets are required.
#'
#' @section GL options:
#' Both option lists accept the following positive finite scalar settings.
#' Orders and capacity controls must be integers. These options are checked
#' against the selected rule; unknown or misspelled names are errors.
#'
#' | Option | Default | Meaning |
#' |---|---|---|
#' | `reference_nodes` | `max(2, floor(nodes/2))` | Distinct GL comparison order, from 2 to 256; may be above or below the fitted order. For observation options, use `observation_nodes` in the default. |
#' | `bound` | 12 | Integration-domain half-width in standard-normal coordinates. The domain remains fixed; a failed omitted-tail check names this option to increase. |
#' | `panel_loglik_tol` | 1e-7 | Per-panel log-likelihood refinement tolerance. |
#' | `panel_score_tol` | 1e-5 | Per-panel first-score refinement tolerance, with coordinate-specific scales. |
#' | `tail_tol` | 1e-9 | Relative omitted probability-tail tolerance. |
#' | `total_loglik_tol` | 1e-5 | Aggregate log-likelihood discrepancy budget across panels. |
#' | `total_score_tol` | 1e-4 | Aggregate first-score discrepancy budget across panels. |
#' | `information_tol` | 2e-3 | Step/order stability tolerance for requested observed information and, when valid, its inverse and SEs. Observation options contribute only for PLN. |
#' | `seconds_per_panel` | 60 | Wall-time limit for one panel's integration work in a preparation or evaluation. |
#' | `seconds_total` | 300 | Wall-time limit across panels in a preparation or evaluation. |
#' | `max_nodes` | 1e8 | Integration-point cap for a panel rule, including tensor-product work where needed. |
#' | `max_cells` | 4096 | Capacity limit for partition refinement on staged routes. |
#' | `batch_nodes` | 8192 outer; 128 observation | Target node batch size for streamed calculations; a complete cell can exceed it. Controls temporary storage, not accuracy. |
#'
#' Score budgets account for coordinate units: coefficient budgets use each
#' regressor's root mean square. A small canonical/reference discrepancy is an
#' empirical interior-error estimate; analytic omitted-tail envelopes are
#' separate. More interior nodes cannot resolve an omitted-tail bound outside
#' the chosen domain. Increase the named resource or domain limit explicitly
#' when justified, then repeat the comparison without relaxing its tolerances.
#' Observation time limits accumulate work across a person's observed rows,
#' including observation-mesh refinement. They are not a fresh allowance for
#' each row. Time guards run between bounded numerical blocks, so an indivisible
#' block can finish after its deadline before the failure is raised. With outer
#' GH/AGH and observation GL, the overall guard is the larger of 300 seconds and
#' the supplied observation total limit. See the staged-quadrature mathematical
#' chapter for the full budget accounting.
#'
#' @section Optimizer methods:
#' | Backend | Accepted methods |
#' |---|---|
#' | `"optim"` | `"BFGS"` (default), `"Nelder-Mead"`, `"CG"`, `"L-BFGS-B"`, `"SANN"` |
#' | `"maxLik"` | `"BFGS"` (default), `"NR"`/`"Newton-Raphson"`, `"BFGSR"`/`"BFGS-R"`, `"BHHH"`, `"CG"`/`"Conjugate-Gradient"`, `"NM"`/`"Nelder-Mead"`, `"SANN"` |
#'
#' Optimizer coordinates are unrestricted coefficients, log SDs and Fisher
#' correlations; supplied starts and reported estimates remain natural-scale.
#' `optim` minimizes the negative log likelihood. The backend status reports
#' that method's termination condition. `SANN` has no convergence test, and
#' a zero-iteration run is reported as not run; neither is labelled converged.
#' Newton methods may require curvature during optimization even with OPG
#' covariance. No method guarantees a global maximum.
#'
#' Without explicit scales, BFGS uses the reciprocal initial panel-score
#' column norms, `1 / sqrt(colSums(score^2))`, for regression coefficients with
#' non-negligible initial scores. Other coordinates use
#' `1 / pmax(1, sqrt(colSums(score^2)))`; all scales are bounded below by
#' `sqrt(.Machine$double.eps)`. This respects coefficient units under rescaled
#' regressors and Linear outcomes. Fixed-GH panel binary selection with free
#' `cor_error` and otherwise default `optim` BFGS settings instead uses
#' `parscale = 1/sqrt(n_panels)` with `reltol = 1e-12`. Explicit `parscale` or
#' `reltol` prevents that automatic pair. The random-effects auxiliaries used
#' by `init = "default"` have fresh BFGS optimizer settings; the joint fit's
#' optimizer options do not control them. Numerical auxiliaries inherit the
#' equation's person-effect integration settings, and PLN also inherits the
#' observation settings. The internal Linear auxiliary uses an exact Gaussian
#' likelihood without quadrature. See [panel_select()] for supplied-start
#' precedence and initialization fallback.
#'
#' The fit's `control` keeps requested settings and `optimizer$control` records
#' effective backend options. `update()` retains the fitted rules and settings
#' unless explicitly replaced; a fresh `panel_control(...)` is a complete
#' control replacement, not a partial patch to the old control.
#'
#' @return A `panelselect_control` list of validated settings, printed by
#'   `print()`.
#' @family fitting functions
#' @seealso [check_quadrature()]
#' @examples
#' panel_control()
#' panel_control(nodes = c(selection = 8, outcome = 12))
#' panel_control(quadrature = "agh", nodes = 15)
#' panel_control(quadrature = "gl", quadrature_options = list(bound = 16))
#' panel_control(optimizer = "optim", method = "L-BFGS-B", factr = 1e6)
#' panel_control(stationarity_warning = NULL)
#' @export
panel_control = function(quadrature = c("gh", "agh", "gl"), nodes = NULL,
                          quadrature_options = list(),
                          observation_quadrature = c("gh", "agh", "gl"),
                          observation_nodes = NULL, observation_options = list(),
                          mode_maxit = 50L, mode_tol = 1e-9, exact_reductions = TRUE,
                          optimizer = c("optim", "maxLik"), method = "BFGS", maxit = 500L,
                          reltol = sqrt(.Machine$double.eps), factr = 1e7, parscale = NULL,
                          trace = 0L, stationarity_warning = 0.1, ...) {
  # `...` stays last so that the named arguments above keep partial matching;
  # a setting given there stops with its replacement.
  earlier = list(...)
  labels = names(earlier)
  if (length(earlier) && (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)))
    stop("panel_control arguments must be uniquely named", call. = FALSE)
  moved = intersect(labels, fit_settings)
  if (length(moved))
    stop(paste(moved, collapse = ", "), if (length(moved) == 1L) " is" else " are",
         " now ", if (length(moved) == 1L) "an argument" else "arguments", " of panel_select()",
         call. = FALSE)
  renamed = intersect(labels, names(renamed_control_arguments))
  if (length(renamed))
    stop(paste0(renamed, " is now ", renamed_control_arguments[renamed], collapse = "; "),
         call. = FALSE)
  removed = intersect(labels, names(removed_control_arguments))
  if (length(removed))
    stop(paste(removed_control_arguments[removed], collapse = "; "), call. = FALSE)
  if (length(labels))
    stop("Unknown panel_control argument(s): ", paste(labels, collapse = ", "),
         "; see ?panel_control", call. = FALSE)
  quadrature = rule_name(quadrature, c("gh", "agh", "gl"), "quadrature")
  observation_quadrature = rule_name(observation_quadrature, c("gh", "agh", "gl"),
    "observation_quadrature")
  optimizer = match.arg(optimizer)
  integer_value = function(value, name, lengths = 1L, minimum = 1L) {
    if (!is.numeric(value) || !length(value) %in% lengths || any(!is.finite(value)) ||
        any(value < minimum | value != floor(value) | value > .Machine$integer.max))
      stop(name, " must contain valid integer values", call. = FALSE)
    as.integer(value)
  }
  if (is.null(nodes)) nodes = if (quadrature == "gl") 20L else 10L
  # A length-two vector orders the selection effect, then the outcome effect, on
  # every route. Names are accepted on input; the stored pair is unnamed.
  if (!is.null(names(nodes))) {
    if (length(nodes) != 2L || !setequal(names(nodes), c("selection", "outcome")))
      stop("Named nodes must be c(selection = , outcome = )", call. = FALSE)
    nodes = unname(nodes[c("selection", "outcome")])
  }
  nodes = integer_value(nodes, "nodes", 1:2)
  quadrature_options = rule_options(quadrature, nodes, quadrature_options, "quadrature_options", 8192L)
  if (is.null(observation_nodes)) observation_nodes = 20L
  observation_nodes = integer_value(observation_nodes, "observation_nodes")
  observation_options = rule_options(observation_quadrature, observation_nodes, observation_options,
    "observation_options", 128L)
  maxit = integer_value(maxit, "maxit", minimum = 0L)
  mode_maxit = integer_value(mode_maxit, "mode_maxit")
  trace = integer_value(trace, "trace", minimum = 0L)
  for (name in c("reltol", "factr", "mode_tol")) {
    value = get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      stop(name, " must be a positive finite number", call. = FALSE)
  }
  if (!is.null(parscale) && (!is.numeric(parscale) || is.complex(parscale) ||
      !length(parscale) || any(!is.finite(parscale) | parscale <= 0)))
    stop("parscale must be NULL or positive finite numeric scales", call. = FALSE)
  if (!is.null(stationarity_warning) && (!is.numeric(stationarity_warning) ||
      is.complex(stationarity_warning) || length(stationarity_warning) != 1L ||
      !is.finite(stationarity_warning) || stationarity_warning <= 0))
    stop("stationarity_warning must be NULL or one positive finite number", call. = FALSE)
  method = optimizer_method(method, optimizer)
  if (!is.logical(exact_reductions) || length(exact_reductions) != 1L || is.na(exact_reductions))
    stop("exact_reductions must be TRUE or FALSE", call. = FALSE)
  structure(list(quadrature = quadrature, nodes = nodes, quadrature_options = quadrature_options,
                 observation_quadrature = observation_quadrature,
                 observation_nodes = observation_nodes, observation_options = observation_options,
                 mode_maxit = mode_maxit, mode_tol = mode_tol, exact_reductions = exact_reductions,
                 optimizer = optimizer, method = method, maxit = maxit, reltol = reltol,
                 factr = factr, parscale = parscale, trace = trace,
                 stationarity_warning = stationarity_warning),
            class = "panelselect_control", reltol_supplied = !missing(reltol))
}

#' @rdname panel_control
#' @export
print.panelselect_control = function(x, ...) {
  x = current_control(x)
  observation = switch(x$observation_quadrature %||% "gh",
    gh = paste("fixed Gauss-Hermite (gh),", x$observation_nodes, "nodes"),
    agh = paste("adaptive Gauss-Hermite (agh),", x$observation_nodes, "nodes"),
    gl = paste("staged Gauss-Legendre (gl), order", x$observation_nodes),
    gk = removed_gk_label)
  options = x$quadrature_options
  outer = switch(x$quadrature,
    gl = paste0("Gauss-Legendre (gl), order ", x$nodes, " per segment (reference order ",
      options$reference_nodes, ", bound ", options$bound, ")"),
    gk = removed_gk_label,
    paste0(if (x$quadrature == "agh") "adaptive" else "fixed", " Gauss-Hermite (", x$quadrature,
      "), ", paste(x$nodes, collapse = " x "), " nodes",
      if (length(x$nodes) == 2L) " (selection x outcome)"))
  cat("PanelSelect control\n")
  cat("Integration: ", outer, if (identical(x$binary_integration, "checked"))
    "; checked binary alias (removed)", if (!isTRUE(x$exact_reductions))
    "; exact reductions off", "\n", sep = "")
  if (identical(x$quadrature, "agh"))
    cat("             mode search: maxit ", x$mode_maxit, ", tolerance ", format(x$mode_tol), "\n", sep = "")
  cat("             PLN observation rule: ", observation, "\n", sep = "")
  cat("Optimizer:   ", x$optimizer, " ", x$method, ", maxit ", x$maxit,
      if (identical(x$method, "L-BFGS-B")) paste0(", factr ", format(x$factr)) else
        paste0(", reltol ", format(x$reltol)),
      if (!is.null(x$parscale)) paste0(", parscale ", paste(format(x$parscale), collapse = " ")),
      if (isTRUE(x$trace > 0)) paste0(", trace ", x$trace), "\n", sep = "")
  cat("Q_OPG warning: above ",
      if (is.null(x$stationarity_warning)) "NULL (warning disabled)" else format(x$stationarity_warning),
      "\n", sep = "")
  # A fitted model's control also stores the fit settings given to panel_select().
  if (length(intersect(fit_settings, names(x))))
    cat("Fit settings: na_action ", x$na_action %||% "-", ", init ", x$init %||% "-", ", vcov ",
        x$vcov %||% "-", "\n", sep = "")
  if (!is.null(x$fixed_geometry)) cat("Frozen adaptive geometry (fixed_geometry, removed)\n")
  if (!is.null(x$integration_mesh))
    cat("Retained staged mesh: stage ", x$integration_mesh$epoch, "\n", sep = "")
  invisible(x)
}

# 2. Input names, removed settings and option validation

# Earlier control names. New calls that use them stop with the replacement;
# stored controls are renamed by current_control().
old_rule_values = c(adaptive = "agh", fixed = "gh", split = "gl")
renamed_control_fields = c(inner_quadrature = "observation_quadrature",
  inner_nodes = "observation_nodes", reduction = "exact_reductions",
  stationarity_tol = "stationarity_warning")
renamed_control_arguments = c(renamed_control_fields, split = "quadrature_options",
  gk = "quadrature_options", inner_split = "observation_options", inner_gk = "observation_options")
renamed_options = c(loglik_tol = "panel_loglik_tol", score_tol = "panel_score_tol",
  aggregate_loglik_tol = "total_loglik_tol", aggregate_score_tol = "total_score_tol",
  max_intervals = "max_cells")

# Removed compatibility settings. New calls stop with the replacement;
# stored controls drop their default values (current_control()).
removed_control_arguments = c(
  binary_integration = 'binary_integration was removed; use quadrature = "gl"',
  checked = 'checked was removed with binary_integration; use quadrature = "gl" and quadrature_options',
  fixed_geometry = 'fixed_geometry was removed; use quadrature = "agh"')

rule_name = function(value, choices, label) {
  if (is.character(value) && length(value) == 1L && value %in% names(old_rule_values))
    stop(label, ' = "', value, '" is now "', old_rule_values[[value]], '"', call. = FALSE)
  if (identical(value, "gk")) stop(label, ' = "gk" was removed; use "gl"', call. = FALSE)
  match.arg(value, choices)
}

# Removed settings: the staged Gauss-Kronrod 21/10 rule, the checked binary
# alias and frozen adaptive geometry. Controls and fits
# saved with them still display (GK with this label) but stop where they would
# integrate, naming the refit.
removed_gk_label = "Gauss-Kronrod 21/10 (gk, removed)"

stop_removed_settings = function(control) {
  if (identical(control$binary_integration, "checked"))
    stop('binary_integration = "checked" was removed; refit with ',
         'update(fit, control = panel_control(quadrature = "gl"))', call. = FALSE)
  if (!is.null(control$fixed_geometry))
    stop("frozen adaptive geometry (fixed_geometry) was removed; refit with ",
         'update(fit, control = panel_control(quadrature = "agh"))', call. = FALSE)
  outer = identical(control$quadrature, "gk")
  inner = identical(control$observation_quadrature, "gk")
  if (!outer && !inner) return(invisible(NULL))
  # The suggested control keeps the other rule when it is not the default.
  rules = c(quadrature = if (outer) "gl" else control$quadrature,
    observation_quadrature = if (inner) "gl" else control$observation_quadrature %||% "gh")
  rules = rules[rules != c("agh", "gh")]
  stop("GK integration was removed; refit with update(fit, control = panel_control(",
       paste0(names(rules), ' = "', rules, '"', collapse = ", "), "))", call. = FALSE)
}

# Optimizer options of the legacy wrappers that no panel_control() argument
# holds: fits keep evaluating, but a refit stops.
stop_legacy_optimizer_options = function(control) {
  options = attr(control, "legacy_optimizer_control", exact = TRUE)
  if (length(options))
    stop("legacy optimizer options (", paste(names(options), collapse = ", "), ") were removed",
         "; refit with update(fit, control = panel_control()), adding the settings to keep",
         call. = FALSE)
}

# Options of the GL rule; the Gauss-Hermite rules take none.
rule_options = function(rule, nodes, options, label, default_batch) {
  if (rule %in% c("agh", "gh")) {
    if (length(options))
      stop(label, " require ", if (label == "quadrature_options") "quadrature" else
        "observation_quadrature", " = 'gl'", call. = FALSE)
    return(list())
  }
  gl_options(nodes, options, default_batch, label)
}

# Option rules below compare full optim names, so abbreviations are expanded.
# maxLik accepts its listed spellings case-insensitively; keep them as given.
optimizer_method = function(method, optimizer) {
  if (!is.character(method) || length(method) != 1L || is.na(method))
    stop("method must be one optimizer method name", call. = FALSE)
  choices = if (optimizer == "optim") c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN") else
    c("NR", "Newton-Raphson", "BFGS", "BFGSR", "BFGS-R", "BHHH", "CG",
      "Conjugate-Gradient", "NM", "Nelder-Mead", "SANN")
  index = if (optimizer == "optim") pmatch(method, choices) else
    match(tolower(method), tolower(choices))
  if (is.na(index))
    stop("method must be one of ", paste(choices, collapse = ", "), " for optimizer = '",
         optimizer, "'", call. = FALSE)
  if (optimizer == "optim") choices[index] else method
}

# Merge a uniquely named list of positive finite scalars into the defaults of
# one rule; `label` names the control argument that holds the list.
merge_positive_options = function(options, defaults, label, rule) {
  if (!is.list(options) || (length(options) && (is.null(names(options)) ||
      anyNA(names(options)) || any(!nzchar(names(options))) || anyDuplicated(names(options)))))
    stop(label, " must be a uniquely named list", call. = FALSE)
  # Staged rules validate their options on every evaluation, so the names are
  # checked with %in% (names are unique here) rather than set operations.
  supplied = names(options)
  old = supplied[supplied %in% names(renamed_options)]
  if (length(old))
    stop(paste0(label, " entry ", old, " is now ", renamed_options[old], collapse = "; "),
         call. = FALSE)
  unknown = supplied[!supplied %in% names(defaults)]
  if (length(unknown))
    stop("Unknown ", label, " entry for the ", rule, " rule: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  for (name in names(options)) {
    value = options[[name]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0)
      stop(label, " ", name, " must be one positive finite real value", call. = FALSE)
    defaults[[name]] = value
  }
  defaults
}

# A selected outer GL rule or PLN observation GL rule excludes posterior APIs.
uses_nonhistory_quadrature = function(control, prepared) {
  uses_gl_quadrature(control) ||
    (identical(prepared$family, "pln") && identical(control$observation_quadrature, "gl"))
}

# GL uses one canonical finite rule; a distinct order only checks it.
gl_options = function(nodes, options = list(), default_batch = 8192L, label = "quadrature_options") {
  if (!is.numeric(nodes) || is.complex(nodes) || length(nodes) != 1L ||
      !is.finite(nodes) || nodes < 3 || nodes > 256 || nodes != floor(nodes))
    stop("GL nodes must be one integer from 3 through 256", call. = FALSE)
  if ("max_bound" %in% names(options))
    stop(label, " entry max_bound belonged to the removed GK rule; GL keeps its bound", call. = FALSE)
  defaults = list(reference_nodes = max(2L, floor(nodes/2)), bound = 12,
    panel_loglik_tol = 1e-7, panel_score_tol = 1e-5, tail_tol = 1e-9,
    total_loglik_tol = 1e-5, total_score_tol = 1e-4,
    information_tol = 2e-3, seconds_per_panel = 60, seconds_total = 300,
    max_nodes = 100000000L, max_cells = 4096L, batch_nodes = default_batch)
  defaults = merge_positive_options(options, defaults, label, "gl")
  if (defaults$reference_nodes != floor(defaults$reference_nodes) ||
      defaults$reference_nodes < 2 || defaults$reference_nodes > 256 || defaults$reference_nodes == nodes)
    stop(label, " reference_nodes must be an integer from 2 through 256 and different from nodes",
         call. = FALSE)
  for (name in c("max_nodes", "max_cells", "batch_nodes")) {
    value = defaults[[name]]
    if (value != floor(value) || value > .Machine$integer.max)
      stop(label, " ", name, " must be a representable positive integer", call. = FALSE)
    defaults[[name]] = as.integer(value)
  }
  defaults$reference_nodes = as.integer(defaults$reference_nodes)
  defaults
}

uses_gl_quadrature = function(control) identical(control$quadrature, "gl")

# 3. Effective controls and normalization of saved controls

# Fit settings belong to panel_select(); a fitted model stores them in its
# control, so revalidation sets them aside.
fit_settings = c("na_action", "init", "vcov")

# The control a fit evaluates: the validated numerical settings with the fit
# settings attached, where the estimation code reads them. The fit settings
# precede a retained mesh, the layout current_control() gives stored controls,
# and the list is built at its final length: serialize() also records a list
# grown in place.
fit_control = function(control = panel_control(), na_action = "fail", init = "default",
                       vcov = "opg") {
  init = validate_initialization(init)
  control = validate_panel_control(control)
  kept = attributes(control)[setdiff(names(attributes(control)), "names")]
  fields = unclass(control)
  mesh = fields$integration_mesh
  fields$integration_mesh = NULL
  answer = c(fields, list(na_action = na_action, init = init, vcov = vcov),
    if (!is.null(mesh)) list(integration_mesh = mesh))
  attributes(answer) = c(list(names = names(answer)), kept)
  answer
}

# Earlier stored controls use rule values adaptive, fixed and split, the
# inner_* fields, one option list per rule (split, gk, inner_split, inner_gk)
# with the earlier option names, reduction and stationarity_tol. They are
# renamed into the current layout and field order, with the fit settings and
# any retained mesh (see current_mesh()) after the numerical settings; current
# controls return at once. A control built by panel_control() that was given
# earlier names afterwards stops: renaming would silently overwrite them.
# Saved controls may also hold the compatibility fields
# fixed_geometry, binary_integration and checked, and legacy fits an attribute
# legacy_optimizer_control: default values (and an empty attribute) are
# dropped; other values stay, so that the control displays them and stops
# where it would use them (stop_removed_settings(),
# stop_legacy_optimizer_options()).
current_control = function(control) {
  lists = c("split", "gk", "inner_split", "inner_gk")
  if (!is.list(control)) return(control)
  # Engine entries call this for every evaluation, so the test for a current
  # control reads fields without copying, dispatch or set operations.
  stored = names(control)
  earlier = stored[stored %in% c(names(renamed_control_fields), lists)]
  old_rule = function(field) {
    value = .subset2(control, field)
    is.character(value) && length(value) == 1L && value %in% names(old_rule_values)
  }
  old_rules = c("quadrature", "observation_quadrature")[c(old_rule("quadrature"),
    old_rule("observation_quadrature"))]
  removed = stored[stored %in% names(removed_control_arguments)]
  options = attr(control, "legacy_optimizer_control", exact = TRUE)
  empty_options = is.list(options) && !length(options)
  parscale = .subset2(control, "parscale")
  old_parameters = any(names(parscale) %in% names(renamed_parameters))
  mesh = .subset2(control, "integration_mesh")
  if (is.list(mesh)) mesh = current_mesh(mesh)
  old_mesh = !identical(mesh, .subset2(control, "integration_mesh"))
  if (!length(earlier) && !length(old_rules) && !length(removed) && !empty_options) {
    # Controls already in the current layout may still need parameter-label
    # translation in a named scale or retained mesh. Keep the fast path
    # and do not rebuild settings or inspect live integration environments.
    if (!old_parameters && !old_mesh) return(control)
    if (old_parameters) names(control$parscale) = current_names(names(parscale))
    if (old_mesh) control$integration_mesh = mesh
    return(control)
  }
  kept = attributes(control)[setdiff(names(attributes(control)),
    c("names", if (empty_options) "legacy_optimizer_control"))]
  answer = unclass(control)
  if (length(earlier) || length(old_rules)) {
    if (any(c(renamed_control_fields, "quadrature_options", "observation_options") %in% names(control)))
      stop("control mixes earlier and current names: ", paste(c(
        if (length(earlier)) paste0(earlier, " is now ", renamed_control_arguments[earlier]),
        if (length(old_rules)) paste0(old_rules, ' = "', unlist(unclass(control)[old_rules]), '" is now "',
          old_rule_values[unlist(unclass(control)[old_rules])], '"')), collapse = "; "), call. = FALSE)
    rule = function(value)
      if (length(value) == 1L && value %in% names(old_rule_values)) old_rule_values[[value]] else value
    option_names = function(value) {
      if (!length(value)) return(list())
      earlier = names(value) %in% names(renamed_options)
      names(value)[earlier] = renamed_options[names(value)[earlier]]
      value
    }
    for (old in intersect(names(renamed_control_fields), names(answer)))
      names(answer)[names(answer) == old] = renamed_control_fields[[old]]
    if (!is.null(answer$quadrature)) answer$quadrature = rule(answer$quadrature)
    if (!is.null(answer$observation_quadrature))
      answer$observation_quadrature = rule(answer$observation_quadrature)
    outer = switch(answer$quadrature %||% "", gl = answer$split, gk = answer$gk, NULL)
    inner = switch(answer$observation_quadrature %||% "", gl = answer$inner_split,
      gk = answer$inner_gk, NULL)
    answer[intersect(lists, names(answer))] = NULL
    answer$quadrature_options = option_names(outer %||% list())
    answer$observation_options = option_names(inner %||% list())
  }
  if (is.null(answer$fixed_geometry)) answer["fixed_geometry"] = NULL
  if (is.null(answer$binary_integration) || identical(answer$binary_integration, "gh"))
    answer["binary_integration"] = NULL
  if (!length(answer$checked)) answer["checked"] = NULL
  if (old_parameters) names(answer$parscale) = current_names(names(parscale))
  if (is.list(mesh)) answer$integration_mesh = mesh
  fields = c(names(panel_control()), fit_settings)
  order = c(intersect(fields, names(answer)), setdiff(names(answer), fields))
  answer = answer[order]
  attributes(answer) = c(list(names = order), kept)
  answer
}

# Earlier saved fits: the control in the earlier vocabulary (see
# current_control()), the component loglik_by_panel, the integration labels
# (rule values adaptive, fixed and split; inner_nodes and inner_quadrature) and
# the person-effect SD names sd_selection and sd_outcome where they index
# parameters: coefficients, covariance matrices, score columns, restrictions,
# the parameter layout, numeric optimizer and integration records, the name map
# stored by removed wrappers and the outcome scale. Labels derived
# from the data (panel IDs, row and variable names) are left alone. Auxiliary
# fits nested in a wrapper's result are read the same way. Other saved fits
# use sd_id_selection, sd_id_outcome and cor_id. Their parameter labels,
# including retained information-check matrices, change without recomputing
# any stored numerical result. Current fits return at once.
# 4. Saved fitted objects and stationarity records

current_fit = function(object) {
  control = current_control(object$control)
  earlier_stationarity = any(c("score_norm", "quadratic") %in%
    names(object$optimizer[["stationarity"]]))
  if (identical(control, object$control) && !"loglik_by_panel" %in% names(object) &&
      !any(names(object$coefficients) %in% names(renamed_parameters)) && !earlier_stationarity)
    return(object)
  named = function(x, parameter_rows = TRUE) {
    if (!is.null(names(x))) names(x) = current_names(names(x))
    if (!is.null(dim(x)) && !is.null(dimnames(x))) {
      axes = seq_along(dimnames(x))
      if (!parameter_rows) axes = axes[axes != 1L]
      for (axis in axes) if (!is.null(dimnames(x)[[axis]]))
        dimnames(x)[[axis]] = current_names(dimnames(x)[[axis]])
    }
    accuracy = attr(x, "checked_accuracy", exact = TRUE)
    if (is.list(accuracy)) attr(x, "checked_accuracy") = records(accuracy)
    x
  }
  # Parameter-indexed numbers and the rule labels of diagnostics, at any depth.
  records = function(x, parameter_rows = TRUE) {
    if (is.data.frame(x)) return(x)
    if (is.list(x)) {
      for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] = records(x[[i]], parameter_rows)
      earlier = names(x) %in% c("inner_nodes", "inner_quadrature")
      names(x)[earlier] = renamed_control_fields[names(x)[earlier]]
      # GL/GK option lists reported by the diagnostics (moving-cut GL: options).
      for (field in intersect(c("options", "inner_options"), names(x))) if (is.list(x[[field]])) {
        earlier = names(x[[field]]) %in% names(renamed_options)
        names(x[[field]])[earlier] = renamed_options[names(x[[field]])[earlier]]
      }
      for (field in intersect(c("method", "observation_quadrature"), names(x))) {
        value = x[[field]]
        if (is.character(value) && length(value) == 1L && value %in% names(old_rule_values))
          x[[field]] = old_rule_values[[value]]
      }
      return(x)
    }
    if (is.numeric(x)) named(x, parameter_rows) else x
  }
  nested = function(x) {
    if (inherits(x, "panelselect")) return(current_fit(x))
    if (is.list(x) && !is.data.frame(x))
      for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] = nested(x[[i]])
    x
  }
  object$control = control
  names(object)[names(object) == "loglik_by_panel"] = "loglik_by_unit"
  for (name in c("coefficients", "fixed")) if (!is.null(object[[name]])) object[[name]] = named(object[[name]])
  if (is.list(object$covariance))
    object$covariance[] = lapply(object$covariance, function(m) if (is.null(m)) m else named(m))
  if (is.matrix(object$score)) colnames(object$score) = current_names(colnames(object$score))
  if (is.data.frame(object$parameter_info))
    object$parameter_info$name = current_names(object$parameter_info$name)
  if (!is.null(object$optimizer)) object$optimizer = records(object$optimizer)
  # Integration matrices index independent units in their first dimension;
  # those row labels are data IDs, even if they spell an old parameter name.
  if (!is.null(object$integration)) object$integration = records(object$integration, parameter_rows = FALSE)
  if (!is.null(object$legacy$parameter_names))
    object$legacy$parameter_names = named(object$legacy$parameter_names)
  if (!is.null(object$prepared$outcome_scale))
    object$prepared$outcome_scale = current_names(object$prepared$outcome_scale)
  for (name in names(object))
    if (!is.null(object[[name]]) && (inherits(object[[name]], "panelselect") || name == "initialization"))
      object[[name]] = nested(object[[name]])
  if (earlier_stationarity) {
    # Only the diagnostic convention changes; retained estimates/information
    # are historical values. Explicit stored cutoffs now use the Q scale.
    optimizer = object$optimizer
    tolerance = if ("stationarity_warning" %in% names(control)) control$stationarity_warning else 0.1
    screen = opg_stationarity(optimizer$gradient, optimizer$covariance$opg, tolerance)
    screen$Q_Hessian = if (isTRUE(optimizer$information_computed))
      hessian_stationarity(optimizer$gradient, optimizer$information,
        optimizer$covariance$observed)$Q_Hessian else NA_real_
    object$optimizer$stationarity = screen
  }
  object
}

# Revalidation of a stored control for a fit, refit or refinement.
validate_panel_control = function(control) {
  control = current_control(control)
  stop_removed_settings(control)
  stop_legacy_optimizer_options(control)
  supplied = attr(control, "reltol_supplied", exact = TRUE) %||% ("reltol" %in% names(control))
  mesh = control$integration_mesh
  arguments = unclass(control)
  arguments[c("integration_mesh", ".rules", fit_settings)] = NULL
  answer = do.call(panel_control, arguments)
  if (!is.null(mesh)) {
    if (!is.list(mesh)) stop("Retained integration_mesh must be a mesh list", call. = FALSE)
    answer$integration_mesh = mesh
  }
  attr(answer, "reltol_supplied") = supplied
  answer
}
