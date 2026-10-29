# Cross-sectional selection with a binary outcome.
#
#   D_i = 1{w_i'alpha + o^s_i + xi_i > 0},  Y*_i = 1{x_i'beta + o^y_i + eps_i > 0},
#   observed when D_i = 1;  (xi_i, eps_i) ~ N2(0, [1 tau; tau 1])
#
#   symbol  parameter    meaning
#   alpha   selection.*  selection coefficients
#   beta    outcome.*    outcome coefficients
#   tau     cor_error    correlation of the shocks
#
# A = a = w'alpha + o^s and B = b = x'beta + o^y. An unselected row contributes
# Phi(-A); a selected row contributes Phi2(A, s*B; s*tau), s = 2y-1. These are
# joint observed-data probabilities, so the selected mass is not divided by
# the selection probability. Missing outcomes are never evaluated.
#
# Read probit_probit_panel for the two row kinds and natural first-score chains.
# It reuses the selected-row partials in ProbitRE_ProbitRE.R; model_assembly.R
# supplies the common output layout. Each row is its own independent unit,
# stored as a one-row panel, with integration dimension q = 0.
# GH controls use integrate_prepared without a person-effect integral; GL uses
# the staged route with no person-effect mesh (quadrature_staged.R). Tau fixed
# at 0 separates the two probit factors; a free tau at 0 retains its derivative.

# Each row is its own unit and there is no person effect: log L_i = log Phi(-a_i)
# (d_i = 0) or log Phi2(a_i, q_i b_i; q_i tau), q_i = 2y_i-1 (d_i = 1), the panel observation
# probability of ProbitRE_ProbitRE.R at A = a, B = b.
probit_probit_panel = function(z, m, out) {
  unselected = node_rows(out, m$data$unselected)
  out = add_rows(out, NULL, local_univariate_probit(m$a[unselected$row], -1, 'a', out$order),
    m$data$unselected, A = list(alpha = 'w'))
  selected = node_rows(out, m$data$selected)
  A = m$a[selected$row]; B = m$b[selected$row]
  if (m$tau_zero) {
    out = add_rows(out, NULL, local_univariate_probit(A, 1, 'a', out$order),
      m$data$selected, A = list(alpha = 'w'))
    return(add_rows(out, NULL, probit_re_rows(B, selected$y, NULL, m, out$order),
      m$data$selected, B = list(beta = 'x')))
  }
  add_rows(out, NULL, probit_re_probit_re_rows(A, B, selected$y, m$tau, out$order),
    m$data$selected, A = list(alpha = 'w'), B = list(beta = 'x'), tau = list(tau = 1))
}

probit_probit_predictions = function(context) {
  response = list(fit = stats::pnorm(context$b), gradient = context$gradient_b*stats::dnorm(context$b))
  if (!context$joint_needed) return(list(response = response))
  joint = prediction_bvn(context$a, context$b, rep(context$tau, length(context$a)), context$gradient_a, context$gradient_b, context$gradient_tau)
  C = exp(joint$log_fit-context$selection$log_fit)
  list(response = response, joint = list(fit = joint$fit, gradient = joint$gradient),
    selected = list(fit = C, gradient = (joint$log_gradient-context$selection$log_gradient)*C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_probit_model = function() model_definition(
  id = 'probit_probit', file = 'Probit_Probit.R', structure = 'cross_sectional',
  family = 'probit', has_selection = TRUE, parameters = 'cor_error',
  panel = probit_probit_panel, predictions = probit_probit_predictions,
  validate_design = validate_probit_probit_design)

validate_probit_probit_design = function(prepared, info) {
  if (!identical(prepared$structure, "cross_sectional") || prepared$family != "probit")
    return(invisible(NULL))
  free = info$name[!info$fixed]
  constant = function(x) length(x) > 0L && all(x == x[1L])
  if (identical(colnames(prepared$w), "(Intercept)") &&
      identical(colnames(prepared$x), "(Intercept)") &&
      all(c("outcome.(Intercept)", "cor_error") %in% free) &&
      constant(prepared$offset_selection) &&
      constant(prepared$offset_outcome[prepared$selected]))
    stop("An intercept-only cross-sectional binary selection model cannot identify ",
         "both the outcome intercept and cor_error; fix one or supply informative ",
         "covariates or known offsets", call. = FALSE)
}
