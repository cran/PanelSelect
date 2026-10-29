#' PanelSelect: panel sample-selection models
#'
#' PanelSelect provides a unified and flexible framework for sample selection 
#' models with continuous, binary, and count outcomes in panel and cross-sectional 
#' data. In particular, it supports a broad family of panel selection models, 
#' which extend the classic Heckman framework to account for persistent 
#' unobserved heterogeneity. The package supports selection at two levels: 
#' the individual level, through correlated random effects, and the observation 
#' level, through correlated occasion-specific shocks. This allows researchers 
#' to distinguish non-random selection driven by persistent individual 
#' characteristics from selection driven by temporary shocks. Alongside model 
#' estimation, PanelSelect supports prediction of population-level outcomes.
#'
#' @section Main functions:
#' - [panel_select()] fits a model from an outcome formula, a selection formula
#'   (or `NULL` for a standalone model), the data, the panel ID column (or
#'   `NULL` for cross-sectional data) and the outcome family; its arguments
#'   also set restrictions, starting values, missing-data handling and the
#'   reported covariance. [panel_control()] sets the integration rules and the
#'   optimizer.
#' - Standard methods ([panelselect-methods]) summarize, compare and update fits;
#'   [predict.panelselect()], [residuals.panelselect()] and [random_effects()]
#'   give predictions and posterior summaries.
#' - [check_quadrature()] compares likelihoods and scores at retained estimates
#'   under higher orders or an independent GL rule; agreement is evidence of
#'   stability, not an integration certificate.
#' - [check_stationarity()] reads stored scores and information, reporting
#'   `Q_OPG` and, with valid positive-definite stored observed information and
#'   its inverse, `Q_Hessian`. The default warning is strictly `Q_OPG > 0.1`.
#'   It does not compute a Hessian or alter backend convergence.
#' - [add_information()] explicitly returns a fit with observed information
#'   stored for reuse. Requested covariance can remain unavailable; optimizer
#'   convergence, stationarity, integration and inference require distinct checks.
#'
#' Use `panel_select()` in place of the PanelSelect 1.x and PanelCount wrappers;
#' `vignette("migration", package = "PanelSelect")` maps their calls and
#' results to `panel_select()`.
#'
#' @section Quick start:
#' ```
#' fit = panel_select(y ~ x, observed ~ .outcome + z, data, id = "person",
#'                    family = "linear")
#' summary(fit)
#' check_stationarity(fit)[c("Q_OPG", "Q_Hessian")]
#' predict(fit, newdata, type = "selected_mean")
#' panel_select(count ~ x, NULL, data, "person", "poisson")  # no selection
#' ```
#'
#' Start with `vignette("vignette", package = "PanelSelect")`; users of earlier
#' versions should read `vignette("migration", package = "PanelSelect")`.
#' GL supports population predictions; history-conditioned predictions and
#' random-effect summaries require a GH/AGH fit, with no prediction-time
#' integration override. History-conditioned prediction requires explicit
#' history; [random_effects()] defaults to the retained fitting rows.
#'
#' @importFrom Rcpp evalCpp
#' @importFrom stats AIC BIC coef confint dnorm fitted formula logLik model.frame
#'   model.matrix na.action nobs pnorm predict qnorm residuals runif rnorm setNames
#'   terms update vcov
#' @importFrom utils head tail
#' @useDynLib PanelSelect, .registration = TRUE
"_PACKAGE"
