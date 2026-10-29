# Cross-sectional selection with a Linear outcome.
#
#   D_i = 1{w_i'alpha + o^s_i + xi_i > 0},  Y*_i = x_i'beta + o^y_i + sigma*eps_i,
#   observed when D_i = 1;  (xi_i, eps_i) ~ N2(0, [1 tau; tau 1])
#
#   symbol  parameter    meaning
#   alpha   selection.*  selection coefficients
#   beta    outcome.*    outcome coefficients
#   sigma   sd_error     outcome SD
#   tau     cor_error    correlation of the shocks
#
# A = a = w'alpha + o^s, B = b = x'beta + o^y and R = (y-B)/sigma. An unselected
# row contributes Phi(-A); a selected row contributes the joint density
#   phi(R)/sigma * Phi((A + tau*R)/sqrt(1-tau^2)).
# Read probit_linear_panel for the row kinds and natural first-score chains.
# It reuses linear_selection_log_contribution (ProbitRE_LinearRE.R) with the
# outcome density retained; model_assembly.R supplies the common output layout.
#
# Each row is an independent unit, stored as a one-row panel, with integration
# dimension q = 0. GH controls use integrate_prepared without quadrature; GL
# uses the staged route without a person-effect mesh (quadrature_staged.R).
# Unselected rows use local_univariate_probit and never evaluate missing outcomes.

# Each row is its own unit: log L_i = log Phi(-a_i) (d_i = 0) or, with
# R_i = (y_i - b_i)/sigma, log phi(R_i) - log sigma + log Phi((a_i + tau R_i)/c_tau).
probit_linear_panel = function(z, m, out) {
  unselected = node_rows(out, m$data$unselected)
  out = add_rows(out, NULL, local_univariate_probit(m$a[unselected$row], -1, 'a', out$order),
    m$data$unselected, A = list(alpha = 'w'))
  selected = node_rows(out, m$data$selected)
  add_rows(out, NULL, linear_selection_log_contribution(m$a[selected$row], m$b[selected$row],
      selected$y, m$error, m$tau, out$order, density = TRUE), m$data$selected,
    A = list(alpha = 'w'), B = list(beta = 'x'), error = list(sigma = 1), tau = list(tau = 1))
}

probit_linear_predictions = function(context) {
  response = list(fit = context$b, gradient = context$gradient_b)
  if (!context$joint_needed) return(list(response = response))
  covariance = context$error*context$tau; gradient_covariance = context$tau*context$gradient_error+
      context$error*context$gradient_tau
  density = stats::dnorm(context$a); mills = prediction_mills(context$a)
  J = context$b*context$selection$fit+covariance*density
  gradient_J = context$gradient_b*context$selection$fit+
      context$selection$gradient*context$b+gradient_covariance*density-
      context$gradient_a*(covariance*context$a*density)
  C = context$b+covariance*mills
  gradient_C = context$gradient_b+gradient_covariance*mills-
      context$gradient_a*(covariance*mills*(context$a+mills))
  list(response = response, joint = list(fit = J, gradient = gradient_J),
       selected = list(fit = C, gradient = gradient_C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_linear_model = function() model_definition(
  id = 'probit_linear', file = 'Probit_Linear.R', structure = 'cross_sectional',
  family = 'gaussian', has_selection = TRUE, parameters = c('sd_error', 'cor_error'),
  panel = probit_linear_panel, predictions = probit_linear_predictions,
  positive_scales = 'sd_error')
