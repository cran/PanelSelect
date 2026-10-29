# Cross-sectional selection with a Poisson-lognormal outcome.
#
#   D_i = 1{w_i'alpha + o^s_i + xi_i > 0},
#   Y*_i | e_i ~ Poisson(exp(x_i'beta + o^y_i + gamma*e_i)), observed when D_i = 1;
#   (xi_i, e_i) ~ N2(0, [1 tau; tau 1])
#
#   symbol  parameter    meaning
#   alpha   selection.*  selection coefficients
#   beta    outcome.*    outcome coefficients
#   gamma   sd_error     observation-level SD
#   tau     cor_error    correlation of xi and e
#
# A = a = w'alpha + o^s and B = b = x'beta + o^y. An unselected row contributes
# Phi(-A); a selected row contributes the joint probability
#   int Phi((A + tau*e)/sqrt(1-tau^2)) Pois(y; exp(B + gamma*e)) phi(e) de.
# Read probit_pln_panel for the row kinds and natural first-score chains. Its
# observation kernel and local partials come from ProbitRE_PLNRE.R; the shared
# observation engine reduces them before model_assembly.R expands the design.
#
# Each row is an independent unit, stored as a one-row panel. There is no
# person effect or outer integral (q = 0), but each selected row needs its own
# observation integral. m$inner uses fixed GH by default, or adaptive GH/staged
# GL as requested (quadrature_inner.R). GL outer controls use the staged route
# without a person-effect mesh. Unselected rows never evaluate missing outcomes.

# The observation kernel is ProbitRE_PLNRE.R's (its selection_normal and
# normal_evaluate hooks are unused here).
probit_pln_inner_kernel = probit_re_pln_re_inner_kernel

# Each row is its own unit: log L_i = log Phi(-a_i) (d_i = 0) or the observation
# integral log int Phi((a_i + tau e)/c_tau) Pois(y_i; e^{b_i + gamma e}) phi(e) de;
# with tau or gamma fixed at 0 the selection factor separates.
probit_pln_panel = function(z, m, out) {
  unselected = node_rows(out, m$data$unselected)
  out = add_rows(out, NULL, local_univariate_probit(m$a[unselected$row], -1, 'a', out$order),
    m$data$unselected, A = list(alpha = 'w'))
  selected = node_rows(out, m$data$selected)
  A = m$a[selected$row]; B = m$b[selected$row]
  log_factorial = m$data$log_factorial[selected$row]
  if (!m$tau_zero && !is.null(m$inner))
    return(add_rows(out, NULL, probit_re_pln_re_rows(A, B, selected$y, log_factorial, m, out$order),
      m$data$selected, A = list(alpha = 'w'), B = list(beta = 'x'), error = list(gamma = 1),
      tau = list(tau = 1)))
  out = add_rows(out, NULL, local_univariate_probit(A, 1, 'a', out$order), m$data$selected,
    A = list(alpha = 'w'))
  add_rows(out, NULL, pln_re_rows(B, selected$y, log_factorial, m, out$order), m$data$selected,
    B = list(beta = 'x'), error = list(gamma = 1))
}

probit_pln_predictions = function(context) {
  logM = context$b+context$error^2/2; gradient_log_M = context$gradient_b+context$error*context$gradient_error
  M = exp(logM); response = list(fit = M, gradient = gradient_log_M*M)
  if (!context$joint_needed) return(list(response = response))
  z = context$a+
      context$error*context$tau; gradient_z = context$gradient_a+
      context$tau*context$gradient_error+
      context$error*context$gradient_tau
  logJ = logM+stats::pnorm(z, log.p = TRUE)
  gradient_log_J = gradient_log_M+gradient_z*prediction_mills(z)
  J = exp(logJ); C = exp(logJ-context$selection$log_fit)
  list(response = response, joint = list(fit = J, gradient = gradient_log_J*J),
    selected = list(fit = C, gradient = (gradient_log_J-context$selection$log_gradient)*C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_pln_model = function() model_definition(
  id = 'probit_pln', file = 'Probit_PLN.R', structure = 'cross_sectional',
  family = 'pln', has_selection = TRUE, parameters = c('sd_error', 'cor_error'),
  panel = probit_pln_panel, predictions = probit_pln_predictions,
  zero_correlations = list(cor_error = 'sd_error'))
