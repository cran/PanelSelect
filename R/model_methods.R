#' Work with fitted PanelSelect models
#'
#' Use the standard R generics below to inspect or reuse a [panel_select()] fit.
#' Reported coefficients and covariance matrices use natural parameter scales.
#'
#' | Task | Functions and help |
#' |---|---|
#' | Extract estimates and uncertainty | [coef()][panelselect-inference], [vcov()][panelselect-inference], [confint()][panelselect-inference] |
#' | Read a fit or coefficient table | [print()][summary.panelselect], [summary()][summary.panelselect] |
#' | Inspect likelihood and compare models | [logLik()][panelselect-comparison], [nobs()][panelselect-comparison], [AIC()][panelselect-comparison], [BIC()][panelselect-comparison] |
#' | Recover formulas, designs and omitted rows | [formula()][panelselect-design], [terms()][panelselect-design], [model.frame()][panelselect-design], [model.matrix()][panelselect-design], [na.action()][panelselect-design] |
#' | Refit with changed settings or data | [update()][update.panelselect] |
#' | Predict or summarize person effects | [predict()][predict.panelselect], [fitted()][predict.panelselect], [residuals()][residuals.panelselect], [random_effects()] |
#' | Check a fit | [check_stationarity()], [check_quadrature()], [add_information()] |
#'
#' See `vignette("vignette", package = "PanelSelect")` for a complete analysis
#' workflow and `vignette("migration", package = "PanelSelect")` for earlier APIs.
#' @name panelselect-methods
#' @seealso [panel_select()], [panel_control()]
NULL

# 1. Fit metadata and argument validation

# A fitted object without the structure field still records its model identity,
# which fixes the structure; missing metadata never selects one.
model_structure = function(object) {
  recorded = object[["structure"]] %||% object$prepared[["structure"]]
  if (!is.null(recorded)) return(recorded)
  model_id = object[["model_id"]] %||% object$prepared[["model_id"]]
  if (is.null(model_id))
    stop("The fitted object has missing metadata: neither a structure nor a model identity", call. = FALSE)
  resolve_model(model_id = model_id)$structure
}

# Methods whose generics pass `...` reject misspelled arguments rather than
# silently using defaults, e.g. summary(fit, type = "observed").
reject_unused = function(dots, method, hint) {
  if (!length(dots)) return(invisible(NULL))
  labels = names(dots) %||% character(length(dots))
  labels[!nzchar(labels)] = "<unnamed>"
  stop(method, "() does not use the argument(s) ", paste(labels, collapse = ", "), "; ", hint,
       call. = FALSE)
}

# 2. Estimates and inference

#' Extract coefficients, covariance matrices and confidence intervals
#'
#' Inspect the natural-scale estimates of a [panel_select()] fit, select a
#' covariance estimator, or obtain parameter confidence intervals. These methods
#' do not refit the model. Observed and sandwich covariance can require an
#' additional information calculation at the fitted estimates.
#'
#' @param object A fitted `panelselect` object from [panel_select()].
#' @param type Covariance estimator for `vcov()`: `"opg"`, `"observed"`, or
#'   `"sandwich"`. Omit it to use the covariance choice stored in the fit.
#'   See the table below for the calculation and when additional work is needed.
#' @param vcov In `confint()`, an estimator name (`"opg"`, `"observed"`, or
#'   `"sandwich"`) or a finite, named natural-scale covariance matrix; the default
#'   is the fit's choice. In `vcov()`, this is an alternative argument name for
#'   `type` and accepts only an estimator name. If both names are supplied to
#'   `vcov()`, their values must agree.
#' @param parm Parameters to include in `confint()`: a character vector of names,
#'   all-positive inclusion indices, all-negative exclusion indices, or a logical
#'   vector with one entry per coefficient. Omit to include every parameter.
#'   Use `names(coef(object))` to find the available names.
#' @param level Confidence level strictly between zero and one; default `0.95`.
#' @param ... `coef()` ignores extra arguments. `vcov()` and `confint()` reject
#'   unused arguments, including misspelled covariance or interval options.
#'
#' @section Choosing a covariance estimator:
#' Let each independent unit contribute a score vector, let `B` be the sum of
#' their outer products, and let `A` be the observed information. Independent
#' units are people in a panel model and rows in a cross-sectional model.
#'
#' | Choice | Calculation | Common use and computation |
#' |---|---|---|
#' | `"opg"` | Inverse of `B` | Default likelihood-based covariance; uses the stored unit scores and needs no Hessian. |
#' | `"observed"` | Inverse of `A` | Curvature-based likelihood inference; computes information if it is not already stored. |
#' | `"sandwich"` | `solve(A) %*% B %*% solve(A)` | Independent-unit sandwich inference; also requires valid observed information. It does not correct an inaccurate likelihood integral or establish model identification. |
#' | Named matrix | Supplied natural covariance | Available to `confint()`, [summary()], and [predict()], for externally calculated uncertainty. See [predict.panelselect()] for matrix validation. |
#'
#' Assign `fit = add_information(fit)` to retain information for repeated
#' observed/sandwich requests. Information that is singular, indefinite, or
#' numerically unavailable is not replaced by a different estimator. Staged GL
#' uses the retained mesh; if information needs a refined mesh, explicitly
#' update/refit with `vcov = "observed"` or `"sandwich"`.
#'
#' @section Parameter scales and intervals:
#' - Regression coefficients are reported on the model's coefficient scale.
#' - SDs are positive natural SDs; intervals are constructed on the log scale
#'   and transformed back. Testing an SD against zero is a boundary problem,
#'   so summaries do not supply a usual zero-test p-value for SDs.
#' - Correlations are reported between -1 and 1; intervals and summary tests
#'   use the Fisher (`atanh`) scale and transform interval limits back.
#' - Fixed parameters remain in `coef()` and have zero covariance rows/columns.
#'   Their inferential entries, including confidence limits, are `NA`.
#'
#' These are Wald intervals. Stored optimizer-coordinate covariance is used when
#' available; supplied matrices and older fits without it use the natural-scale
#' delta method. Numerical convergence alone does not establish interval coverage.
#'
#' @return `coef()` returns a named numeric vector. `vcov()` returns a square
#'   named natural-scale covariance matrix. `confint()` returns a two-column
#'   matrix of lower and upper limits for the requested parameters.
#' @examples
#' dat = data.frame(person = rep(1:20, each = 3), y = rep(c(0, 1, 3, 1, 2, 4), 10))
#' fit = panel_select(y ~ 1, NULL, dat, "person", "poisson",
#'                    fixed = c(sd_re_outcome = 0))
#' coef(fit)
#' vcov(fit)
#' confint(fit, parm = "outcome.(Intercept)")
#' fit = add_information(fit)
#' vcov(fit, type = "observed")
#' @rdname panelselect-inference
#' @aliases panelselect-inference
#' @seealso [summary.panelselect()], [add_information()], [check_stationarity()]
#' @export
coef.panelselect = function(object, ...) current_fit(object)$coefficients

#' @rdname panelselect-inference
#' @export
vcov.panelselect = function(object, type, vcov, ...) {
  object = current_fit(object)
  reject_unused(list(...), "vcov",
    'the estimator is selected with type = or vcov = ("opg", "observed" or "sandwich")')
  # vcov = is the name summary(), confint(), predict() and panel_select() use.
  if (!missing(type) && !missing(vcov) && !identical(type, vcov))
    stop("type and vcov both name the covariance estimator; give one of them", call. = FALSE)
  estimator = if (!missing(type)) type else if (!missing(vcov)) vcov else object$control$vcov
  covariance_pair(object, estimator)$natural
}

#' @rdname panelselect-inference
#' @export
confint.panelselect = function(object, parm, level = 0.95, vcov = object$control$vcov, ...) {
  object = current_fit(object)
  reject_unused(list(...), "confint", "it takes parm =, level = and vcov =")
  table = coefficient_table(object, vcov, level)
  parameters = rownames(table)
  if (missing(parm)) parm = parameters else if (is.character(parm)) {
    earlier = intersect(parm, names(renamed_parameters))
    if (length(earlier))
      stop(paste0(earlier, " is now ", renamed_parameters[earlier], collapse = "; "), " (in parm)",
           call. = FALSE)
    if (any(!parm %in% parameters))
      stop("Unknown parameter in parm: ", paste(setdiff(parm, parameters), collapse = ", "), call. = FALSE)
  } else if (is.logical(parm)) {
    if (length(parm) != length(parameters) || anyNA(parm))
      stop("A logical parm must give one TRUE or FALSE value per parameter", call. = FALSE)
    parm = parameters[parm]
  } else if (is.numeric(parm) && !anyNA(parm) && all(parm == floor(parm)) &&
             all(abs(parm) <= length(parameters)) && (all(parm > 0) || all(parm < 0))) {
    parm = parameters[parm]
  } else stop("parm must contain parameter names, positive or negative indices, or one ",
              "logical value per parameter", call. = FALSE)
  result = table[parm, c("lci", "uci"), drop = FALSE]
  colnames(result) = paste0(format(100 * c((1-level)/2, (1+level)/2), trim = TRUE), " %")
  result
}

# 3. Likelihood and model comparison

#' Inspect likelihood, sample size and information criteria
#'
#' Extract the fitted joint log likelihood or compare models fitted to the same
#' sample and response definition. Within-person dependence is included in each
#' panel's integrated likelihood contribution; observations within a person are
#' not treated as independent likelihood units.
#'
#' @param object A fitted `panelselect` model. For a multiple-model AIC/BIC table,
#'   this is the first model to compare.
#' @param ... Additional fitted models in `AIC()` and `BIC()`. They must supply
#'   compatible `logLik()` and `nobs()` methods and represent comparable
#'   likelihoods. Extra arguments are ignored by `logLik()` and `nobs()`.
#' @param k Penalty per free parameter in `AIC()`; default `2` gives ordinary
#'   AIC. Supply one finite number. Fixed parameters do not contribute to the
#'   parameter count.
#' @param unit Sample-size convention for `BIC()`: `"panels"` (default) uses
#'   independent people for panel models and independent rows for cross-sectional
#'   models; `"observations"` uses retained rows for every model.
#'
#' @section Calculations and sample counts:
#' If `ell` is the fitted log likelihood and `p` is the number of free parameters,
#' `AIC = -2 * ell + k * p` and `BIC = -2 * ell + log(n) * p`.
#'
#' | Quantity or option | Count |
#' |---|---|
#' | `nobs(fit)` | Retained rows, including observed nonparticipation rows in a selection model. |
#' | `attr(logLik(fit), "nobs")` | Independent people for a panel model; retained rows for a cross-sectional model. |
#' | `attr(logLik(fit), "n_rows")` | Retained rows. |
#' | `attr(logLik(fit), "df")` | Estimated coefficients, SDs and correlations; fixed parameters excluded. |
#' | `BIC(fit, unit = "panels")` | Uses the independent-unit `nobs` attribute of `logLik()`. |
#' | `BIC(fit, unit = "observations")` | Uses `nobs(fit)` instead. |
#'
#' @section Comparing models responsibly:
#' - Use the same retained sample and response definition. Equal sample counts
#'   do not establish sample identity, and criteria from different likelihood
#'   definitions need not be comparable.
#' - The existing guard warns for different row counts, and for default BIC
#'   also different independent-unit counts. It does not compare row identities.
#' - Check integration sensitivity and optimization before interpreting small
#'   criterion differences. A higher-order fit need not be a more accurate fit.
#' - Use likelihood-ratio tests only with the appropriate nesting, identified
#'   parameters and boundary conditions; see [panel_select()].
#'
#' @return `logLik()` returns a `logLik` scalar with `df`, `nobs` and `n_rows`
#'   attributes. `nobs()` returns the retained row count. AIC/BIC return one
#'   number for one model or a data frame with `df` and the criterion for several.
#' @examples
#' dat = data.frame(person = rep(1:20, each = 3), y = rep(c(0, 1, 3, 1, 2, 4), 10))
#' fit = panel_select(y ~ 1, NULL, dat, "person", "poisson",
#'                    fixed = c(sd_re_outcome = 0))
#' logLik(fit)
#' c(rows = nobs(fit), people = attr(logLik(fit), "nobs"))
#' AIC(fit)
#' BIC(fit)
#' BIC(fit, unit = "observations")
#' @rdname panelselect-comparison
#' @aliases panelselect-comparison
#' @seealso [panel_select()], [update.panelselect()], [check_quadrature()]
#' @export
logLik.panelselect = function(object, ...) {
  units = if (model_structure(object) == "cross_sectional") object$prepared$n_rows else
    object$prepared$n_panels
  structure(object$loglik, df = sum(!object$parameter_info$fixed),
            nobs = units, n_rows = object$prepared$n_rows,
            class = "logLik")
}

#' @rdname panelselect-comparison
#' @export
nobs.panelselect = function(object, ...) object$prepared$n_rows

model_criterion = function(object, dots, call, bic = FALSE, k = 2, unit = "panels") {
  models = c(list(object), dots)
  likelihoods = lapply(models, stats::logLik)
  df = vapply(likelihoods, attr, numeric(1), which = "df")
  # Count checks catch some sample mismatches; equal counts do not establish
  # that the same rows or response definition were used.
  units = vapply(likelihoods, function(x) as.numeric(attr(x, "nobs") %||% NA), numeric(1))
  rows = vapply(models, stats::nobs, numeric(1))
  if (length(dots) && (length(unique(rows)) > 1L ||
                       (bic && unit == "panels" && length(unique(units)) > 1L)))
    warning("models are not all fitted to the same sample (rows or independent units differ)",
            call. = FALSE)
  penalty = if (bic) log(if (unit == "panels") units else rows) else k
  value = -2 * vapply(likelihoods, as.numeric, numeric(1)) + penalty * df
  if (!length(dots)) return(unname(value))
  labels = vapply(c(list(call$object), as.list(call$...)), deparse1, character(1))
  result = data.frame(df = df, value = value, row.names = make.unique(labels))
  names(result)[2L] = if (bic) "BIC" else "AIC"
  result
}

#' @rdname panelselect-comparison
#' @export
AIC.panelselect = function(object, ..., k = 2) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k))
    stop("k must be one finite number", call. = FALSE)
  model_criterion(object, list(...), match.call(expand.dots = FALSE), k = k)
}

#' @rdname panelselect-comparison
#' @export
BIC.panelselect = function(object, ..., unit = c("panels", "observations")) {
  model_criterion(object, list(...), match.call(expand.dots = FALSE), bic = TRUE,
                  unit = match.arg(unit))
}

# 4. Formulas and designs

#' Recover fitted formulas, designs and omitted rows
#'
#' Inspect either equation of a fitted model using its retained data, factor
#' levels and contrasts. The methods describe the fitted design rather than
#' constructing a design for new prediction data.
#'
#' @param x,object A fitted `panelselect` model. `formula()` and `terms()` call
#'   this argument `x`; `model.matrix()` and `na.action()` call it `object`.
#' @param formula A fitted `panelselect` model, despite the generic argument name,
#'   when calling `model.frame()`.
#' @param equation Equation to extract: `"outcome"` (default) or `"selection"`.
#'   A standalone model has no selection equation: `formula()` returns `NULL`,
#'   while `terms()`, `model.frame()` and `model.matrix()` stop for that request.
#' @param ... Extra generic arguments; these methods do not use them. In
#'   particular, they do not accept new data or a replacement missing-data policy.
#'
#' @section Returned designs:
#' - `formula()` returns the fitted equation formula, with `.outcome` selection
#'   shorthand expanded to the regressors used by the fit.
#' - `terms()` returns its terms object, including offset information.
#' - `model.frame()` rebuilds the retained rows in input order, using the fitted
#'   factor levels and allowing structural missing outcome entries.
#' - `model.matrix()` constructs the corresponding coefficient design with the
#'   fitted contrasts. Offsets belong to the formula/model frame, not additional
#'   estimated-coefficient columns.
#' - `na.action()` identifies rows omitted or excluded for required missing data;
#'   it returns `NULL` when none were dropped. Rows outside `subset` are separate.
#'
#' Outcome values/regressors unused on unselected rows can remain missing in
#' the returned frame. The likelihood's internal zero placeholders are not
#' imputed outcomes. Categories represented only among unselected rows can
#' produce an unidentified outcome coefficient and a rank-deficiency error;
#' levels are not silently dropped as a model-selection step. Prepare an
#' estimable outcome design and check new-data factor support explicitly.
#'
#' @return A formula, terms object, data frame, numeric design matrix, or
#'   `omit`/`exclude` row-index object, according to the generic.
#' @examples
#' dat = data.frame(person = rep(1:20, each = 3), y = rep(c(0, 1, 3, 1, 2, 4), 10))
#' fit = panel_select(y ~ 1, NULL, dat, "person", "poisson",
#'                    fixed = c(sd_re_outcome = 0))
#' formula(fit)
#' head(model.frame(fit))
#' head(model.matrix(fit))
#' na.action(fit)
#' @rdname panelselect-design
#' @aliases panelselect-design
#' @seealso [panel_select()], [predict.panelselect()], [update.panelselect()]
#' @export
formula.panelselect = function(x, equation = c("outcome", "selection"), ...) {
  equation = match.arg(equation)
  x$prepared$design[[equation]]$formula
}

#' @rdname panelselect-design
#' @export
model.frame.panelselect = function(formula, equation = c("outcome", "selection"), ...) {
  equation = match.arg(equation)
  design = formula$prepared$design[[equation]]
  if (is.null(design)) stop("This model has no selection equation", call. = FALSE)
  frame = stats::model.frame(design$terms,
                              strip_contrasts(formula$data[formula$prepared$keep, , drop = FALSE]),
                              na.action = stats::na.pass, xlev = design$xlevels,
                              drop.unused.levels = FALSE)
  attr(frame, "na.action") = formula$prepared$na.action
  frame
}

#' @rdname panelselect-design
#' @export
terms.panelselect = function(x, equation = c("outcome", "selection"), ...) {
  equation = match.arg(equation)
  design = x$prepared$design[[equation]]
  if (is.null(design)) stop("This model has no selection equation", call. = FALSE)
  design$terms
}

#' @rdname panelselect-design
#' @export
na.action.panelselect = function(object, ...) object$prepared$na.action

#' @rdname panelselect-design
#' @export
model.matrix.panelselect = function(object, equation = c("outcome", "selection"), ...) {
  equation = match.arg(equation)
  design = object$prepared$design[[equation]]
  frame = model.frame.panelselect(object, equation)
  stats::model.matrix(design$terms, frame, contrasts.arg = design$contrasts)
}

# 5. Updated fitting calls

#' Refit a model with changed formulas, data or settings
#'
#' Refit a [panel_select()] model using retained data and fitted coefficients as
#' warm starts. Settings and restrictions persist unless explicitly overridden.
#' Use `start = NULL` to run an initializer afresh instead of warm-starting.
#'
#' @param object A fitted `panelselect` model. Its retained data and settings
#'   supply defaults, so the original data-frame variable need not still exist.
#' @param formula. Optional updated outcome formula following
#'   [stats::update.formula()], for example `. ~ . + z` or `. ~ . - x`.
#'   Omit to keep the fitted outcome formula.
#' @param ... Named overrides of [panel_select()] arguments: `outcome`,
#'   `selection`, `data`, `id`, `family`, `start`, `fixed`, `control`, `subset`,
#'   `na_action`, `init`, or `vcov`. Values are evaluated in the calling
#'   environment; `subset` is evaluated in the data. See the options below.
#' @param evaluate Logical; default `TRUE` fits the updated model. `FALSE`
#'   returns a self-contained fitting call after argument checks and design
#'   preparation; rank, identification and initialization checks run when the
#'   returned call is evaluated. This constructs the call without
#'   running initialization or optimization. The returned call can embed data.
#'
#' @section Common updates:
#' | Change | How to request it |
#' |---|---|
#' | Add/remove outcome terms | `update(fit, . ~ . + z)`; a `.outcome` selection shorthand follows the updated outcome. |
#' | Change selection terms | `update(fit, selection = . ~ . + z)` or supply a complete selection formula. |
#' | Override some warm starts | `start = c(cor_re = 0.2)`; other surviving coefficients remain warm starts. |
#' | Try a fresh initializer | `start = NULL, init = "pooled"` (or another supported `init`); changing `init` alone normally leaves complete warm starts in place. |
#' | Change restrictions | Supply a complete replacement `fixed` vector; `fixed = numeric()` removes explicit restrictions. Structural zero-scale restrictions still apply. |
#' | Change one numerical control | `control = modifyList(fit$control, list(maxit = 1000))` preserves the remaining controls. |
#' | Replace all numerical controls | `control = panel_control(...)` uses current defaults for omitted controls, including GH for omitted `quadrature`. |
#' | Change default inference | `vcov = "observed"`, `"opg"`, or `"sandwich"`; this is a fitting argument, not a control-list option. |
#' | Use new data | Supply `data`; if the fit used `subset`, also supply its new value or `subset = NULL`. |
#'
#' @section Starting values and changed models:
#' Explicit starts must be finite named numeric vectors with unique, valid
#' parameter names; matrices, unnamed vectors and duplicate names are rejected
#' before they can overwrite warm starts. Named partial starts remain supported.
#' Parameters absent from a changed design/family leave inherited starts and
#' restrictions. Newly freed SDs previously fixed at zero get automatic starts.
#' An empty start vector adds no overrides; `start = NULL` requests automatic
#' starts for all free parameters. Explicit start/fixed conflicts and nonpositive
#' free SD starts are rejected when the fitting call is evaluated.
#'
#' The stored `na_action`, `init`, `vcov`, integration rule and node orders
#' persist unless overridden. Staged meshes are rebuilt for each update. No
#' method is silently changed to improve convergence or integration accuracy.
#'
#' @section Formula shorthand and saved fits:
#' In a dotted selection update, `.outcome` is one parenthesized term, so
#' `. ~ . - x` cannot remove just one regressor inside that shorthand; write
#' the selection formula in full for that change. Earlier native 2.0 fits are
#' read using current names, but removed integration/optimizer settings stop
#' with guidance for an explicit replacement. See the migration vignette.
#' The former initializer string `"ProbitRE"` can remain in a displayed saved
#' fit, but it is not a supported update/refit setting. Supply `init = "default"`
#' (or another current initializer) explicitly, even with `evaluate = FALSE`.
#' `update(fit, init = "default")` keeps fitted warm starts;
#' `update(fit, init = "default", start = NULL)` reruns the automatic
#' initialization described in [panel_select()].
#' Earlier person-effect parameter names in a saved fit are read as
#' `sd_re_selection`, `sd_re_outcome` and `cor_re`; explicit new `start` and
#' `fixed` inputs must use these current names. The `id` data-column argument
#' is unchanged.
#' Previously stored Gaussian information is not repaired merely by printing or
#' calling [add_information()] again; update/refit explicitly for the corrected
#' outcome-scale differentiation steps.
#'
#' @return A newly fitted `panelselect` model, or a fitting call when
#'   `evaluate = FALSE`. The supplied fit is unchanged.
#' @examples
#' dat = data.frame(person = rep(1:20, each = 3), y = rep(c(0, 1, 3, 1, 2, 4), 10))
#' fit = panel_select(y ~ 1, NULL, dat, "person", "poisson",
#'                    fixed = c(sd_re_outcome = 0))
#' observed = update(fit, vcov = "observed")
#' fresh = update(fit, start = NULL, init = "pooled")
#' call = update(fit, control = modifyList(fit$control, list(maxit = 1000)),
#'               evaluate = FALSE)
#' coef(fresh)
#' @seealso [panel_select()], [panel_control()], [panelselect-methods]
#' @export
update.panelselect = function(object, formula., ..., evaluate = TRUE) {
  object = current_fit(object)
  # The selection keeps its `.outcome` shorthand, so an updated outcome
  # re-expands it; `formula.` and dotted overrides follow update.formula().
  # The fit settings stored with the control become panel_select() arguments,
  # so overriding one does not conflict with the inherited control.
  control = object$control
  arguments = list(outcome = stats::formula(object),
    selection = object$prepared$design$selection$shorthand %||%
      stats::formula(object, equation = "selection"),
    data = object$data,
    id = if (model_structure(object) == "cross_sectional") NULL else object$prepared[["id"]],
    family = object$prepared$family,
    start = object$coefficients, fixed = object$fixed, control = control,
    subset = object$prepared$subset,
    # Controls stored before init existed keep the pooled starts they were fitted with.
    na_action = control$na_action %||% "fail", init = control$init %||% "pooled",
    vcov = control$vcov %||% "opg")
  if (is.list(control)) arguments$control[fit_settings] = NULL
  if (!missing(formula.))
    arguments$outcome = stats::update.formula(arguments$outcome, formula.)
  extras = match.call(expand.dots = FALSE)$...
  if (length(extras) && (is.null(names(extras)) || any(names(extras) == "") ||
                        anyDuplicated(names(extras)) ||
                        any(!names(extras) %in% names(arguments))))
    stop("Update overrides must name panel_select arguments: ",
         paste(names(arguments), collapse = ", "), call. = FALSE)
  if ("data" %in% names(extras) && !is.null(arguments$subset) && !"subset" %in% names(extras))
    stop("This fit used subset rows; supply subset = for the new data, or subset = NULL",
         call. = FALSE)
  for (name in setdiff(names(extras), "subset"))
    arguments[name] = list(eval(extras[[name]], parent.frame()))
  # A subset is evaluated in the (possibly new) data, as panel_select() does.
  if ("subset" %in% names(extras))
    arguments["subset"] = list(eval(extras$subset, as.data.frame(arguments$data), parent.frame()))
  for (name in intersect(c("outcome", "selection"), names(extras))) {
    value = arguments[[name]]
    if (inherits(value, "formula") && "." %in% all.vars(value)) {
      previous = if (name == "outcome") stats::formula(object) else
        object$prepared$design$selection$shorthand %||% stats::formula(object, equation = "selection")
      if (is.null(previous)) stop("This fit has no ", name, " equation to update", call. = FALSE)
      arguments[[name]] = stats::update.formula(previous, value)
    }
  }
  # Retained meshes reproduce saved evaluations; refits prepare their own mesh.
  if (is.list(arguments$control)) arguments$control$integration_mesh = arguments$control$.rules = NULL
  if (model_structure(object) != "cross_sectional" && is.null(arguments$id) &&
      !"id" %in% names(extras))
    stop("A saved panel model needs an explicit id to update; missing metadata does not select a cross-sectional model", call. = FALSE)
  validate_panel_control(arguments$control)
  # match.arg() would read NULL as its first choice: a NULL family or setting
  # must not quietly select the Linear outcome or another default.
  for (name in c("family", "na_action", "init", "vcov")) if (is.null(arguments[[name]]))
    stop(name, " is required in update(): omit it to keep the fitted value", call. = FALSE)
  arguments$init = validate_initialization(arguments$init)
  match.arg(arguments$vcov, eval(formals(panel_select)$vcov))
  family = family_name(arguments$family)
  prepared = prepare_data(arguments$outcome, arguments$selection, arguments$data,
    arguments$id, family, match.arg(arguments$na_action, c("fail", "omit", "exclude")),
    subset = subset_rows(arguments$subset, as.data.frame(arguments$data)))
  info = parameter_layout(prepared)
  if (!"fixed" %in% names(extras))
    arguments$fixed = arguments$fixed[names(arguments$fixed) %in% info$name]
  fixed = resolve_fixed(arguments$fixed, info, prepared)
  # Warm starts are the inherited coefficients in the new layout; supplied
  # starts replace single entries, and start = NULL requests automatic starts.
  supplied = if ("start" %in% names(extras)) arguments$start
  if (!is.null(supplied))
    validate_named_parameters(supplied, info, "start", allow_zero = TRUE)
  if (!"start" %in% names(extras) || !is.null(supplied)) {
    start = object$coefficients[names(object$coefficients) %in% info$name]
    start[intersect(names(start), names(fixed))] = fixed[intersect(names(start), names(fixed))]
    zero_free_sd = names(start) %in% info$name[info$type == "sd"] &
      !names(start) %in% names(fixed) & start == 0
    start = start[!zero_free_sd]
    start[names(supplied)] = supplied
    arguments["start"] = list(start)
  }
  if (!evaluate) return(as.call(c(list(quote(PanelSelect::panel_select)), arguments)))
  fit = do.call("panel_select", arguments, envir = environment(panel_select))
  fit$call = update_record(object, extras, arguments, !missing(formula.))
  fit
}

# A readable record of an update or refit: the parent call in the order of
# the panel_select() arguments, with the overrides as written. Evaluation used
# the resolved values; warm starts and embedded data frames from older records
# are not repeated.
update_record = function(object, extras, arguments, outcome_updated) {
  call = object$call
  if (is.call(call)) {
    # Earlier saved calls may name the removed formula alias, omit the
    # selection of a standalone model, the id of a cross-sectional one or the
    # default family, and give fit settings inside panel_control().
    if ("formula" %in% names(call)) names(call)[names(call) == "formula"] = "outcome"
    if (isFALSE(object$prepared$has_selection) && !"selection" %in% names(call))
      call["selection"] = list(NULL)
    if (model_structure(object) == "cross_sectional" && !"id" %in% names(call)) call["id"] = list(NULL)
    if (!"family" %in% names(call) && !is.null(object$prepared$family)) call$family = object$prepared$family
    control = call$control
    if (is.call(control) && deparse(control[[1L]]) %in% c("panel_control", "PanelSelect::panel_control"))
      for (name in intersect(fit_settings, names(control))) {
        if (!name %in% names(call)) call[[name]] = control[[name]]
        control[[name]] = NULL
        call$control = control
      }
    call = tryCatch(match.call(panel_select, call), error = function(e) NULL)
  }
  if (!is.call(call)) call = quote(PanelSelect::panel_select())
  for (i in seq_along(call)[-1L]) if (is.data.frame(call[[i]])) call[[i]] = as.name("<retained data>")
  if (outcome_updated) call$outcome = arguments$outcome
  for (name in names(extras)) call[name] = list(if (name %in% c("outcome", "selection") &&
    inherits(arguments[[name]], "formula")) arguments[[name]] else extras[[name]])
  call
}
