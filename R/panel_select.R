#' Fit panel and cross-sectional sample-selection models
#'
#' @description
#' `panel_select()` fits panel and cross-sectional sample-selection models
#' with probit selection and linear, binary, or count outcomes. It also
#' supports standalone random-effects panel models for binary and count
#' outcomes.
#'
#' Panel selection models allow dependence between selection and outcomes
#' through correlated individual-level random effects, correlated
#' observation-level errors, or both. These models are described in
#' Bailey and Peng (2025) and Peng and Van den Bulte (2024).
#' 
#' Estimation uses maximum likelihood with analytical gradients and numerical
#' integration of random effects. The mathematical engine
#' is deeply optimized using C++. The table below summarizes the available models.
#' 
#' @param outcome Outcome formula, such as `y ~ x1 + x2`, with
#'   optional `offset()` terms. Binary responses must be zero/one, counts
#'   nonnegative integers, and linear responses finite numeric values on the
#'   rows where the outcome is observed.
#' @param selection Formula for the selection model, such
#'   as `observed ~ x1 + x2 + z` or `observed ~ .outcome + z`, where `.outcome`
#'   stands for all outcome regressors and `z` represents an optional regressor
#'   affecting selection but not outcome. Set to `NULL` for a standalone probit, Poisson,
#'   or Poisson-lognormal random-effects model without selection. `.outcome`
#'   copies the outcome regressor terms, but not its offsets or intercept choice.
#' @param data A data frame containing the variables used by the formulas,
#'   along with a column indicating panel ID for panel models. Outcome observations
#'   on unselected rows are unused and may be missing. Selection regressors are required
#'   on every retained row; see Formulas and missing data.
#' @param id Name of the panel ID column for a panel model. Set `id = NULL`
#'   to fit a cross-sectional sample-selection model. Supplying an ID column
#'   always fits a panel model, even if each ID occurs only once. Standalone
#'   random-effects models (`selection = NULL`) also require an ID column.
#' @param family Outcome model: `"linear"` for linear regression with normal
#'   errors (aliases `"continuous"` and `"gaussian"`), `"binary"` for probit
#'   regression (alias `"probit"`), or `"count"` for Poisson-lognormal
#'   regression (alias `"pln"`). Panel models include person-level random
#'   effects. The Poisson-lognormal model also includes a normal observation-level
#'   effect on the log rate, allowing additional overdispersion and,
#'   in selection models, correlation with the selection error (`cor_error`).
#'   In panel data, the `"poisson"` option becomes available and fits a Poisson model
#'   with person-level random effects, but no observation-level random effects.
#'   This option still permits overdispersion through person-level heterogeneity
#'   and allows selection-outcome
#'   dependence only through correlated person effects (`cor_re`).
#' @param start `NULL` for default starts, or a finite, uniquely named numeric
#'   vector of natural-scale values. Partial vectors replace only those starts;
#'   a complete free-parameter vector bypasses default initialization. Free
#'   SDs must start positive and correlations strictly between -1 and 1. See
#'   Models for parameter names.
#' @param fixed Named parameter values, on their natural scales, to hold fixed
#'   during estimation. For example, `fixed = c(cor_error = 0)` sets the
#'   correlation between the selection and outcome errors within each observation
#'   to zero, which can substantially speed up estimation. Values supplied in
#'   both `start` and `fixed` must agree. See Restrictions with `fixed`.
#' @param control Numerical integration and optimization settings from
#'   [panel_control()]. By default, random effects are integrated using fixed
#'   Gauss-Hermite quadrature with 10 nodes per person effect and 20 nodes for
#'   Poisson-lognormal observation effects. The default optimizer is `optim`
#'   with the BFGS method and a maximum of 500 iterations. For binary panel
#'   outcomes, Gauss-Legendre quadrature is recommended when the outcome
#'   random-effect SD (`sd_re_outcome`) is large (e.g., 5 or higher). It can
#'   improve accuracy, often at greater computational cost. For example,
#'   `control = panel_control(quadrature = "gl", nodes = 15)` uses 15 nodes per
#'   segment with Gauss-Legendre quadrature.
#' @param subset Optional, subset of rows to use for estimation, specified as a
#'   logical vector, row indices, or row names. Expressions are evaluated in
#'   `data`, so column names can be used directly (e.g., `age >= 18`). `NA`
#'   entries in a logical subset are handled according to `na_action`.
#' @param na_action Handling of missing required values: `"fail"` (default)
#'   stops and identifies rows; `"omit"` drops them; `"exclude"` drops them and
#'   restores `NA` positions in fitting-row predictions and residuals. This
#'   does not drop structurally missing unselected outcomes.
#' @param init Initialization method for free parameters not supplied in `start`.
#'   `"default"` (alias `"re"`) uses estimates from separate random-effects
#'   selection and outcome models for panel selection, and from pooled cross-sectional
#'   regressions for standalone and cross-sectional models. `"pooled"` (alias
#'   `"cross_sectional"`) always uses pooled cross-sectional regressions.
#'   `"zero"` sets regression coefficients and correlations to zero,
#'   and standard deviations to one. `"unif"` and `"norm"` draw initial values
#'   from Uniform(0, 1) and Normal(0, 1), respectively. Standard deviations use
#'   `exp(draw)` and correlations use `tanh(draw)`. Use `set.seed()` for
#'   reproducible random starts. See Starting values for details.
#' @param vcov Covariance estimator used for standard errors and inference.
#'   `"opg"` (default) uses the inverse of the summed independent-unit score
#'   outer products and avoids a post-fit Hessian calculation. `"observed"`
#'   uses the inverse observed information (the negative Hessian of the log
#'   likelihood). `"sandwich"` combines the observed information with the
#'   score outer products to provide asymptotically robust standard errors,
#'   clustered by panel ID for panel models. Selecting `"observed"` or
#'   `"sandwich"` requests a numerical Hessian calculation after optimization,
#'   which adds computation time. See Inference for assumptions and details
#'   on computing the Hessian after fitting.
#'
#' @section Models:
#' | Model | `family` | `selection` | `id` | Structural parameters |
#' |---|---|---|---|---|
#' | `ProbitRE_LinearRE` | `"linear"` | Yes | Yes | `sd_re_selection`, `sd_re_outcome`, `sd_error`, `cor_re`, `cor_error` |
#' | `ProbitRE_ProbitRE` | `"probit"` | Yes | Yes | `sd_re_selection`, `sd_re_outcome`, `cor_re`, `cor_error` |
#' | `ProbitRE_PoissonRE` | `"poisson"` | Yes | Yes | `sd_re_selection`, `sd_re_outcome`, `cor_re` |
#' | `ProbitRE_PLNRE` | `"pln"` | Yes | Yes | `sd_re_selection`, `sd_re_outcome`, `sd_error`, `cor_re`, `cor_error` |
#' | `Probit_Linear` | `"linear"` | Yes | No | `sd_error`, `cor_error` |
#' | `Probit_Probit` | `"probit"` | Yes | No | `cor_error` |
#' | `Probit_PLN` | `"pln"` | Yes | No | `sd_error`, `cor_error` |
#' | `ProbitRE` | `"probit"` | No | Yes | `sd_re_outcome` |
#' | `PoissonRE` | `"poisson"` | No | Yes | `sd_re_outcome` |
#' | `PLNRE` | `"pln"` | No | Yes | `sd_re_outcome`, `sd_error` |
#'
#' In the `selection` and `id` columns, Yes means that the argument is required, while No means it must be `NULL`.
#'
#' Regression coefficients use `selection.<column>` and `outcome.<column>`.
#' `sd_re_selection` and `sd_re_outcome` are individual-effect SDs, with
#' correlation `cor_re`. `sd_error` is the residual SD for Linear outcomes or
#' the observation-level random effect SD for PLN outcomes. `cor_error` measures
#' selection-outcome error correlation within an observation.
#' `summary()` describes the structural parameters for the fitted model.
#'
#' @section Formulas and missing data:
#' In selection models, indicators and selection regressors are required on
#' every retained row.
#' Outcomes and outcome regressors are required only on selected rows; unused
#' entries may be `NA`. With `selection = NULL`, outcomes are required on every
#' retained row, and omitting missing outcomes does not correct for selection.
#'
#' Use `.outcome` to copy outcome regressors into the selection formula;
#' `observed ~ .outcome + z` adds an exclusion restriction that can aid
#' identification. Ordinary `.` shorthand is unsupported. Formulas support
#' factors, interactions, transformations and `offset()`.
#'
#' `subset` is applied before missing-data checks. `na_action` handles missing
#' required values: `"fail"` stops, `"omit"` drops rows, and `"exclude"` also
#' restores `NA` positions in fitting-row predictions and residuals.
#' See [panelselect-design] for design and factor-support details.
#'
#'
#' @section Restrictions with `fixed`:
#' `fixed` holds named parameters, including regression coefficients, at
#' natural-scale values. Use the parameter names listed under Models.
#'
#' - `fixed = c(cor_re = 0)` removes dependence between individual effects.
#' - `fixed = c(cor_error = 0)` removes same-observation error dependence and
#'   can substantially speed up estimation. Use
#'   `fixed = c(cor_re = 0, cor_error = 0)` to combine both restrictions when
#'   both correlations are present.
#' - Setting an SD to zero removes that random effect and automatically fixes
#'   any unidentified associated correlation to zero. In PLN models,
#'   `fixed = c(sd_error = 0)` removes observation-level heterogeneity,
#'   reducing the outcome to Poisson.
#'
#' Correlations must remain strictly between -1 and 1; Linear `sd_error` must
#' be positive. Fixed parameters remain in the estimates, have no inferential
#' entries, and do not count toward AIC/BIC.
#'
#' @section Starting values:
#' `init` supplies values not specified in `start`. The default uses separate
#' random-effects fits for applicable panel selection equations; standalone and
#' cross-sectional models use pooled starts. Failed, nonconvergent, or invalid
#' auxiliary fits warn and retain pooled starts for the affected equation.
#' Use `init = "pooled"` to bypass auxiliary fitting; see `init` for other choices.
#'
#' Named partial starts override automatic values; a complete free-parameter
#' vector bypasses initialization. Starts and fixed values use natural scales
#' and must agree where both are supplied.
#'
#' `update()` normally uses fitted coefficients as starts. To initialize afresh,
#' use `update(fit, start = NULL, init = "pooled")`; changing `init` alone normally
#' preserves fitted starts. See [update.panelselect()].
#'
#' @section Inference:
#' Use [summary()][summary.panelselect], [vcov()][panelselect-inference] and
#' [confint()][panelselect-inference] for standard errors and intervals.
#' The default `vcov = "opg"` uses unit scores and avoids a Hessian calculation.
#' `"observed"` uses observed information; `"sandwich"` combines it with unit
#' scores. Both require an additional information calculation.
#'
#' OPG and observed-information inference assume a correctly specified
#' likelihood. Sandwich inference is asymptotic and requires independent units
#' and consistent coefficient estimates; it does not correct coefficient bias
#' or integration error. Independent units are panels, or rows for
#' cross-sectional models.
#'
#' `fit = add_information(fit)` stores observed information without refitting
#' or changing the default covariance. Unavailable requested covariance is
#' never replaced automatically. See [panelselect-inference] for formulas and
#' interval construction, and [add_information()] for repeated requests.
#'
#' @section Integration and numerical checks:
#' [panel_control()] describes GH, AGH and GL settings. GL orders count nodes
#' per integration segment, so equal GH/AGH and GL orders are not directly
#' comparable. PLN's observation-level rule is set separately by
#' `observation_quadrature`; changing the person-effect rule leaves it unchanged.
#'
#' Use [check_stationarity()] to inspect optimizer convergence and `Q_OPG`;
#' the package warns when `Q_OPG > 0.1` by default. Stationarity does not
#' establish integration accuracy. [check_quadrature()] compares likelihoods
#' and scores with parameters fixed. Refit explicitly with `update()` or
#' `refit = TRUE` to assess changes in estimates.
#' Agreement supports numerical stability but does not guarantee accuracy
#' or a global maximum.
#'
#' GL supports population predictions, but history-conditioned predictions
#' and [random_effects()] are unavailable when GL is used at either integration
#' level. See [predict.panelselect()] for supported predictions.
#'
#' @return An object of class `panelselect`. Use
#' [summary()][summary.panelselect] to view the estimation results,
#' including parameter estimates and standard errors.
#' The object contains these main components:
#'
#' - `coefficients`: named regression coefficients, SDs and correlations on
#'   natural scales, including fixed values. Extract them with
#'   [coef()][panelselect-inference].
#' - `covariance`: stored parameter covariance matrices. Use
#'   [vcov()][panelselect-inference] and [summary()][summary.panelselect] for
#'   uncertainty estimates with the chosen estimator.
#' - `loglik`: total log likelihood at the returned estimates; extract with
#'   [logLik()][panelselect-comparison].
#' - `optimizer`: convergence and stationarity information; inspect with
#'   [check_stationarity()].
#' - `integration`: available numerical integration diagnostics. Assess
#'   sensitivity to the integration rule and node counts with [check_quadrature()].
#' - `control` and `fixed`: fitting settings and parameter restrictions.
#'   Modify them through [update()][update.panelselect].
#'
#' See [panelselect-methods] for the full set of fitted-model methods.
#' @references
#' Bailey, M., & Peng, J. (2025). A Random Effects Model of Non-Ignorable Nonresponse
#' in Panel Survey Data. Available at SSRN \doi{10.2139/ssrn.5475626}
#'
#' Peng, J., & Van den Bulte, C. (2024). Participation vs. Effectiveness in Sponsored
#' Tweet Campaigns: A Quality-Quantity Conundrum. *Management Science*, 70(11):7961-7983.
#' Available at \doi{10.1287/mnsc.2019.01897}
#' @family fitting functions
#' @seealso [panelselect-methods], [predict.panelselect()],
#'   `vignette("vignette", package = "PanelSelect")`
#' @examples
#' # 1. Fit a panel selection model
#' set.seed(1)
#' dat = data.frame(person = rep(1:60, each = 4), x = rnorm(240), z = rnorm(240))
#' u = rnorm(60); v = 0.5 * u + rnorm(60, sd = 0.8)  # correlated person effects
#' dat$observed = as.integer(0.3 + 0.5 * dat$z + 0.7 * u[dat$person] + rnorm(240) > 0)
#' dat$y = ifelse(dat$observed == 1, 1 + 0.5 * dat$x + v[dat$person] + rnorm(240), NA)
#'
#' # z enters the selection equation but is excluded from the outcome.
#' fit = panel_select(y ~ x, observed ~ .outcome + z, dat, id = "person",
#'                    family = "linear")
#' summary(fit)
#'
#' # 2. Check and refine the fit
#' check_stationarity(fit)[c("optimizer_converged", "Q_OPG")]
#' # Compare calculations with more nodes, holding the estimates fixed.
#' check_quadrature(fit, nodes = 20)
#' \donttest{
#' # The initial fit uses default GH; refit with 20 nodes.
#' fit20 = update(fit, control = panel_control(nodes = 20))
#' cbind(GH10 = coef(fit), GH20 = coef(fit20))
#' }
#'
#' # 3. Predict population outcome means
#' targets = data.frame(x = c(-1, 0, 1), z = 0)
#' # Mean potential outcomes before selection, averaging over person effects.
#' predict(fit, targets, type = "response",
#'         conditioning = "population", interval = "confidence")
#'
#' \donttest{
#' # 4. Set both selection-outcome correlations to zero
#' independent = update(fit, fixed = c(cor_re = 0, cor_error = 0))
#' coef(independent)
#'
#' # 5. Fit a cross-sectional selection model
#' # Use one observation per independent person.
#' cross_dat = dat[!duplicated(dat$person), ]
#' cross = panel_select(y ~ x, observed ~ .outcome + z, data = cross_dat,
#'                      id = NULL, family = "linear")
#' coef(cross)
#'
#' # 6. Fit a standalone Poisson model with exposure
#' # Counts are observed for everyone, so there is no selection equation.
#' dat$exposure = runif(240, 0.5, 2)
#' dat$count = rpois(240, dat$exposure * exp(0.2 + 0.3 * dat$x + 0.4 * u[dat$person]))
#' counts = panel_select(count ~ x + offset(log(exposure)), selection = NULL,
#'                       data = dat, id = "person", family = "poisson")
#' coef(counts)
#' }
#' @export
panel_select = function(outcome, selection, data, id, family, start = NULL, fixed = NULL,
                         control = panel_control(), subset = NULL,
                         na_action = c("fail", "omit", "exclude"),
                         init = c("default", "re", "pooled", "cross_sectional", "zero", "unif", "norm"),
                         vcov = c("opg", "observed", "sandwich")) {
  call = match.call()
  if (missing(outcome)) stop("outcome is required: the two-sided outcome formula", call. = FALSE)
  if (missing(selection))
    stop("selection is required: give the selection formula (indicator ~ regressors), or ",
         "selection = NULL for a standalone probit, Poisson or PLN random-effects model",
         call. = FALSE)
  # The earlier order was (outcome, data, id, family, selection), so a call
  # in that order puts the data, or with data = named the id, where the
  # selection equation now goes.
  if (!is.null(selection) && !inherits(selection, "formula"))
    stop("panel_select() takes (outcome, selection, data, id, family): the second argument ",
         "is the selection formula or NULL, and data comes third (received ",
         if (is.data.frame(selection)) "a data frame" else paste("an object of class", class(selection)[1L]),
         ")", call. = FALSE)
  if (missing(data)) stop("data is required", call. = FALSE)
  if (missing(id))
    stop("id is required: name the panel ID column, or write id = NULL for independent ",
         "cross-sectional observations", call. = FALSE)
  if (missing(family) || is.null(family)) stop("family is required: ", family_choices, call. = FALSE)
  begin = proc.time()[["elapsed"]]
  # 1. Validate the public arguments and numerical controls.
  family = family_name(family)
  settings = list(na_action = match.arg(na_action),
    init = validate_initialization(if (missing(init)) "default" else init), vcov = match.arg(vcov))
  if (!is.list(control)) stop("control must be a control list", call. = FALSE)
  # A fitted or saved control can also store fit settings; the function
  # arguments determine the settings used.
  carried = intersect(fit_settings, names(control))
  differ = carried[!vapply(carried, function(key) {
    value = control[[key]]
    if (key == "init")
      value = tryCatch(validate_initialization(value), error = function(e) value)
    identical(value, settings[[key]])
  }, TRUE)]
  if (length(differ))
    warning("control carries ", paste0(differ, " = \"", unlist(unclass(control)[differ]), "\"",
      collapse = ", "), ", but these are panel_select() arguments; the fit uses ",
      paste0(differ, " = \"", unlist(settings[differ]), "\"", collapse = ", "), call. = FALSE)
  control = fit_control(control, settings$na_action, settings$init, settings$vcov)
  # Stop before data preparation and auxiliary fits when maxLik is missing.
  if (control$optimizer == "maxLik") maxlik_namespace()
  # 2. Prepare the retained sample, restrictions and natural-scale starts.
  data = as.data.frame(data)
  subset = subset_rows(eval(substitute(subset), data, parent.frame()), data)
  prepared = prepare_data(outcome, selection, data, id, family, control$na_action,
    subset = subset)
  info = parameter_layout(prepared)
  fixed = resolve_fixed(fixed, info, prepared)
  info$fixed = info$name %in% names(fixed)
  validate_design_rank(prepared, info, fixed)
  theta = initial_parameters(prepared, info, start, fixed, control)
  # 3. Optimize the requested finite integration objective.
  if (uses_staged_integration(prepared, control, fixed)) {
    staged = fit_staged_model(theta, info, prepared, control, fixed)
    control = staged$control; evaluator = staged$evaluator
    result = staged$result; final = staged$final; options = staged$options
    optimize_time = staged$optimize_time; rejections = staged$rejections
    information_computed = staged$information_computed
    information_error = staged$information_error; information = staged$information
    evaluations = staged$evaluations
    theta = final$theta
  } else {
    # Preserve the established GH and analytical moving-cut GL fit path.
    eta = parameter_to_optimizer(theta, info)
    evaluator = fit_evaluator(theta, info, prepared, control, fixed)
    initial = evaluator$evaluate(eta, trial = FALSE)
    optimize_begin = proc.time()[["elapsed"]]
    options = optimizer_options(control, prepared, info, initial$optimizer_score)
    result = optimize_model(eta, evaluator, control, options)
    optimize_time = proc.time()[["elapsed"]] - optimize_begin
    final = evaluator$evaluate(result$par, trial = FALSE)
    rejections = if (is.function(evaluator$rejections)) evaluator$rejections() else NULL
    theta = final$theta
    information_computed = control$vcov != "opg"
    information_error = NULL
    information = if (information_computed) information_matrix(result$par, evaluator,
      on_error = function(message) information_error <<- message) else
      matrix(NA_real_, length(eta), length(eta), dimnames = list(names(eta), names(eta)))
    evaluations = if (is.function(evaluator$evaluations)) evaluator$evaluations()
  }
  # 4. Assemble inference and separate convergence from stationarity.
  covariance = covariance_matrices(theta, info, final$score, information)
  status = optimizer_status(result$convergence, control, options, length(result$par))
  converged = status %in% c("converged", "all_fixed")
  gradient = colSums(final$optimizer_score)
  stationarity = opg_stationarity(gradient, covariance$optimizer$opg, control$stationarity_warning)
  stationarity$Q_Hessian = if (information_computed)
    hessian_stationarity(gradient, information, covariance$optimizer$observed)$Q_Hessian else NA_real_
  if (status == "not_converged") warning("Optimizer did not report convergence: ", result$message %||% result$convergence,
                         call. = FALSE)
  if (status == "converged" && stationarity$warning) {
    note = boundary_note(theta, info, resolve_model(prepared$family, prepared$has_selection,
      prepared$structure, prepared$model_id))
    warning("Optimizer stopped with Q_OPG = ", signif(stationarity$Q_OPG, 4),
      " above ", stationarity$tolerance, if (!is.null(note)) paste0("; ", note),
      "; inspect check_stationarity() and integration diagnostics", call. = FALSE)
  }
  if (any(!is.finite(covariance$natural[[control$vcov]][!info$fixed, !info$fixed, drop = FALSE])))
    warning("The requested covariance estimate is unavailable; inspect the fit diagnostics", call. = FALSE)
  # 5. Retain the fitted state and label independent-unit contributions.
  # Panel log likelihoods and scores are named by panel ID; cross-sectional
  # contributions by the row names of their observations.
  units = if (prepared$structure == "panel") as.character(prepared$id_values) else
    rownames(data)[prepared$original_index]
  score = final$score; rownames(score) = units
  structure(list(model_id = prepared$model_id, structure = prepared$structure,
                 coefficients = theta, covariance = covariance$natural,
                 loglik = sum(final$loglik), loglik_by_unit = stats::setNames(final$loglik, units),
                 score = score,
                 optimizer = list(par = result$par,
                                  control = options,
                                  gradient = gradient,
                                  information = information,
                                  information_computed = information_computed,
                                  information_error = information_error,
                                  covariance = covariance$optimizer,
                                  counts = result$counts, code = result$convergence,
                                  message = result$message, converged = converged, status = status,
                                  stationarity = stationarity, rejections = rejections,
                                  time = optimize_time, evaluations = evaluations),
                 integration = staged_panel_summaries(final$diagnostics), prepared = prepared,
                 control = control, fixed = fixed, parameter_info = info,
                 data = data, call = call, time = proc.time()[["elapsed"]] - begin),
            class = "panelselect")
}

family_choices = paste0("\"linear\" (aliases \"continuous\" and \"gaussian\"), ",
  "\"probit\" (alias \"binary\"), \"pln\" (alias \"count\") or \"poisson\"")

# Names and aliases accept unique prefixes and map to one internal family.
# Keep "gaussian" internally so existing model dispatch and saved fits agree.
family_name = function(family) {
  if (!is.character(family) || length(family) != 1L)
    stop("family must be one of ", family_choices, "; a family object such as binomial() is not used",
         call. = FALSE)
  choices = c(linear = "gaussian", continuous = "gaussian", gaussian = "gaussian",
    probit = "probit", binary = "probit", pln = "pln", count = "pln", poisson = "poisson")
  index = pmatch(family, names(choices))
  if (is.na(index))
    stop('family = "', family, '" is not available; use ', family_choices, call. = FALSE)
  unname(choices[index])
}
