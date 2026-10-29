# Panel probit selection with a Poisson outcome.
#
#   D_it = 1{w_it'alpha + o^s_it + delta*u_i + xi_it > 0}
#   Y*_it | v_i ~ Poisson(exp(x_it'beta + o^y_it + lambda*v_i)), observed when D_it = 1
#   (u_i, v_i) ~ N2(0, [1 rho; rho 1]); given (u, v) the count and xi are independent
#
#   symbol  parameter        meaning
#   alpha   selection.*      selection coefficients
#   beta    outcome.*        outcome coefficients
#   delta   sd_re_selection  person-level selection SD
#   lambda  sd_re_outcome    person-level outcome SD
#   rho     cor_re           correlation of the person effects
#
# The model point m holds a_t and b_t. With both person effects active, normal
# nodes z give u = z1, v = rho*z1 + sqrt(1 - rho^2)*z2, so A_t = a_t + delta*u and
# B_t = b_t + lambda*v. Fixed GH written as one loop:
#   for each node k: r_ik = sum_{d=0} log Phi(-A_tk)
#                         + sum_{d=1} [log Phi(A_tk) + y_t B_tk - exp(B_tk) - log y_t!]
#   log L_i = log sum_k w_k exp(r_ik),  score_i = sum_k pi_ik g_i(z_k)  (pi: node posterior)
# Here g_i is the natural score at a fixed node. Read probit_re_poisson_re_panel
# for the two observation factors and their parameter chains; the Poisson
# partials are in poisson_re_rows (PoissonRE.R). Missing outcomes enter only
# log Phi(-A), and model_assembly.R sums rows and assembles mixed derivatives.
#
# Routes (likelihood.R dispatches reductions before this model's entry):
#   fixed GH, both effects       integrate_selection_fixed: Phi(s_t A_t) reused across v | u
#   adaptive GH, history         probit_re_poisson_re_panel via integrate_person_effects
#   rho fixed 0                  independent ProbitRE and PoissonRE factors (likelihood.R)
#   panels without outcomes      ProbitRE selection component (likelihood.R)
#   GL (staged meshes)           probit_re_poisson_re_panel on fixed cells (quadrature_staged.R)
# GH factorization and the nested fixed rule require exact_reductions; the
# nested rule also requires both person effects active and no history output.
# Fixed-zero rho permits factorization; a free rho evaluated at zero keeps its
# coupled score. History retains the unreduced latent representation.

# r_i(z) = sum_{d=0} log Phi(-A_t) + sum_{d=1} [log Phi(A_t) + log Pois(y_t; exp(B_t))]:
# given (u, v) the count and the selection shock are independent.
# A_t = a_t + delta*u, B_t = b_t + lambda*v, u = z1, v = rho*z1 + sqrt(1-rho^2)*z2.
# g_alpha = sum w l_A, g_delta = u sum l_A, g_beta = sum x l_B, g_lambda = v sum l_B,
# g_rho = lambda (dv/drho) sum l_B.
probit_re_poisson_re_panel = function(z, m, out) {
  e = person_effects(m, z, out)
  dA = list(alpha = 'w', delta = e$u)
  dB = list(beta = 'x', lambda = e$v, rho = m$lambda*e$v_rho)
  unselected = node_rows(out, m$data$unselected)
  A = m$a[unselected$row]+m$delta*e$u[unselected$unit]
  out = add_rows(out, e, local_univariate_probit(A, -1, 'a', out$order),
    m$data$unselected, A = dA)
  selected = node_rows(out, m$data$selected)
  A = m$a[selected$row]+m$delta*e$u[selected$unit]
  B = m$b[selected$row]+m$lambda*e$v[selected$unit]
  out = add_rows(out, e, local_univariate_probit(A, 1, 'a', out$order), m$data$selected, A = dA)
  add_rows(out, e, poisson_re_rows(B, selected$y, m$data$log_factorial[selected$row], m,
    out$order), m$data$selected, B = dB)
}

# Likelihood. Phi(s_t A_t) depends on u only (there is no tau), so fixed GH over
# both effects conditions on u and integrates v | u over the outcome nodes
# (integrate_selection_fixed: the same two-dimensional rule, with the selection
# factor reused across outcome nodes). Otherwise integrate_person_effects
# integrates the panel kernel (adaptive GH, history rules).
probit_re_poisson_re_likelihood = function(theta, prepared, control, fixed, q, score, posterior) {
  m = make_model_context(theta, prepared, control, fixed)
  if (q == 2L && control$quadrature != 'agh' && !posterior &&
      isTRUE(control$exact_reductions) && m$tau_zero)
    return(integrate_selection_fixed(theta, prepared, m, control$.rules$outer[[3L]], score))
  integrate_person_effects(theta, prepared, control, fixed, q, score, posterior, m = m)
}

# Person-effect tilting shifts the selection index by delta*lambda*rho.
probit_re_poisson_re_predictions = function(context) {
  mean = count_prediction_mean(context)
  if (!context$joint_needed) return(mean)
  z = context$a+if (context$conditional) 0 else context$delta*context$lambda*context$rho/context$selection_scale
  gradient_z = context$gradient_a+if (context$conditional) context$zero else
    (context$lambda*context$rho*context$gradient_delta+
      context$delta*context$rho*context$gradient_lambda+
      context$delta*context$lambda*context$gradient_rho)/context$selection_scale -
      context$gradient_delta*(context$delta*context$lambda*context$rho*context$delta/context$selection_scale^3)
  logJ = mean$log_mean+stats::pnorm(z, log.p = TRUE)
  gradient_log_J = mean$log_gradient+gradient_z*prediction_mills(z)
  J = exp(logJ); C = exp(logJ-context$selection$log_fit)
  list(response = mean$response, joint = list(fit = J, gradient = gradient_log_J*J),
    selected = list(fit = C, gradient = (gradient_log_J-context$selection$log_gradient)*C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_re_poisson_re_model = function() model_definition(
  id = 'probit_re_poisson_re', file = 'ProbitRE_PoissonRE.R', structure = 'panel',
  family = 'poisson', has_selection = TRUE, parameters = c('sd_re_selection', 'sd_re_outcome', 'cor_re'),
  likelihood = probit_re_poisson_re_likelihood, panel = probit_re_poisson_re_panel,
  outcome_rows = poisson_re_rows, predictions = probit_re_poisson_re_predictions,
  zero_correlations = list(cor_re = c('sd_re_selection', 'sd_re_outcome')))
