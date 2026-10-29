# Panel probit selection with a binary outcome.
#
#   D_it  = 1{w_it'alpha + o^s_it + delta*u_i + xi_it > 0}
#   Y*_it = 1{x_it'beta + o^y_it + lambda*v_i + eps_it > 0}, observed when D_it = 1
#   (u_i, v_i) ~ N2(0, [1 rho; rho 1]),  (xi_it, eps_it) ~ N2(0, [1 tau; tau 1])
#
#   symbol  parameter        meaning
#   alpha   selection.*      selection coefficients
#   beta    outcome.*        outcome coefficients
#   delta   sd_re_selection  person-level selection SD
#   lambda  sd_re_outcome    person-level outcome SD
#   rho     cor_re           correlation of the person effects
#   tau     cor_error        correlation of the observation shocks
#
# The model point m (make_model_context) holds a_t = w_t'alpha + o^s_t and
# b_t = x_t'beta + o^y_t. With both person effects active, independent normal
# nodes z give u = z1 and v = rho*z1 + sqrt(1 - rho^2)*z2, so A_t = a_t + delta*u and
# B_t = b_t + lambda*v. Fixed GH written as one loop (the engine evaluates a
# batch of nodes at once):
#   for each node k: r_ik = sum_{d=0} log Phi(-A_tk) + sum_{d=1} log Phi2(A_tk, q_t B_tk; q_t tau)
#   log L_i = log sum_k w_k exp(r_ik),  score_i = sum_k pi_ik g_i(z_k)  (pi: node posterior)
# q_t = 2y_t - 1 is the response sign; g_i is the natural score at a fixed node.
# Read probit_re_probit_re_rows for the selected joint mass and partials, then
# probit_re_probit_re_panel for selected/unselected sums and parameter chains.
# Missing outcomes enter only log Phi(-A); no outcome partial is evaluated.
#
# Routes (likelihood.R dispatches reductions before this model's entry):
#   fixed/adaptive GH, history   probit_re_probit_re_panel via integrate_person_effects
#   fixed GH, tau fixed 0        integrate_selection_fixed: Phi(s_t A_t) reused across v | u
#   rho = tau = 0 fixed          independent ProbitRE factors (likelihood.R)
#   panels without outcomes      ProbitRE selection component (likelihood.R)
#   GL, tau fixed 0 (moving cuts) split_physical_panel, split_single_outcome_panel
#                                (quadrature_split.R section 3)
#   other GL (staged meshes)     probit_re_probit_re_panel on fixed cells (quadrature_staged.R)
# GH factorization and nested shortcuts require exact_reductions; the nested
# fixed rule also requires both person effects active and no history output.
#
# model_assembly.R applies the declared natural first-score chains; its index
# maps supply latent and higher mixed derivatives. Fixed-zero restrictions
# permit factorization; free correlations evaluated at zero retain their scores.
# History bypasses likelihood reductions to retain the latent prediction target.

# Selected row, q = 2y - 1: log Phi2(A, qB; q tau) and its partials in A, B, tau.
probit_re_probit_re_rows = function(A, B, y, tau, order) {
  sign = 2*y-1
  # The joint mass is Phi_2(A, sign*B; sign*tau), without conditioning on D=1.
  signed_B = sign*B; r = sign*tau
  variance = 1-r*r; c_tau = sqrt(variance)
  logp = log_bvn_value(A, signed_B, r)
  if (!order) return(local_block(logp, c('a', 'b', 'tau'), 2L, order))
  d_A = exp(stats::dnorm(A, log = TRUE) +
    stats::pnorm((signed_B-r*A)/c_tau, log.p = TRUE) - logp)
  d_signed_B = exp(stats::dnorm(signed_B, log = TRUE) +
    stats::pnorm((A-r*signed_B)/c_tau, log.p = TRUE) - logp)
  d_r = exp(-log(2*pi)-log(variance)/2 -
    (A*A-2*r*A*signed_B+signed_B*signed_B)/(2*variance) - logp)
  first = unname(cbind(d_A, d_signed_B, d_r))
  signs = unname(cbind(rep(1, length(A)), sign, sign))
  out = local_block(logp, c('a', 'b', 'tau'), 2L, order, first*signs)
  if (order >= 2L) out = binary_log_higher(out, A, signed_B, r, first, signs, order)
  out
}

# r_i(z) = sum_{d=0} log Phi(-A_t) + sum_{d=1} log Phi2(A_t, q_t B_t; q_t tau),
# A_t = a_t + delta*u, B_t = b_t + lambda*v, u = z1, v = rho*z1 + sqrt(1-rho^2)*z2.
# g_alpha = sum w l_A, g_delta = u sum l_A, g_beta = sum x l_B, g_lambda = v sum l_B,
# g_rho = lambda (dv/drho) sum l_B, g_tau = sum l_tau.
probit_re_probit_re_panel = function(z, m, out) {
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
  if (m$tau_zero) {
    # Fixed tau = 0: Phi2(A, qB; 0) = Phi(A) Phi(qB).
    out = add_rows(out, e, local_univariate_probit(A, 1, 'a', out$order),
      m$data$selected, A = dA)
    return(add_rows(out, e, probit_re_rows(B, selected$y, NULL, m, out$order),
      m$data$selected, B = dB))
  }
  add_rows(out, e, probit_re_probit_re_rows(A, B, selected$y, m$tau, out$order),
    m$data$selected, A = dA, B = dB, tau = list(tau = 1))
}

# Likelihood. With tau fixed at 0 and fixed GH over both effects, Phi(s_t A_t)
# depends on u only, so the rule conditions on u and integrates v | u over the
# outcome nodes (integrate_selection_fixed: the same two-dimensional rule, with
# the selection factor reused across outcome nodes). Otherwise
# integrate_person_effects integrates the panel kernel (fixed or adaptive GH,
# history rules).
probit_re_probit_re_likelihood = function(theta, prepared, control, fixed, q, score, posterior) {
  m = make_model_context(theta, prepared, control, fixed)
  if (q == 2L && control$quadrature != 'agh' && !posterior &&
      isTRUE(control$exact_reductions) && m$tau_zero)
    return(integrate_selection_fixed(theta, prepared, m, control$.rules$outer[[3L]], score))
  integrate_person_effects(theta, prepared, control, fixed, q, score, posterior, m = m)
}

probit_re_probit_re_predictions = function(context) {
  mean = probit_re_predictions(context)
  if (!context$joint_needed) return(mean)
  scale = mean$scale
  correlation = if (context$conditional) context$tau else (context$delta*context$lambda*context$rho+
      context$tau)/(context$selection_scale*scale)
  gradient_correlation = if (context$conditional) context$gradient_tau else
    (context$lambda*context$rho*context$gradient_delta+
      context$delta*context$rho*context$gradient_lambda+
      context$delta*context$lambda*context$gradient_rho+
      context$gradient_tau)/
      (context$selection_scale*scale) - correlation*(context$gradient_delta*context$delta/context$selection_scale^2+
      context$gradient_lambda*context$lambda/scale^2)
  joint = prediction_bvn(context$a, mean$index, rep(correlation, length(context$a)), context$gradient_a, mean$index_gradient, gradient_correlation)
  C = exp(joint$log_fit-context$selection$log_fit)
  list(response = mean$response, joint = list(fit = joint$fit, gradient = joint$gradient),
       selected = list(fit = C, gradient = (joint$log_gradient-context$selection$log_gradient)*C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_re_probit_re_model = function() model_definition(
  id = 'probit_re_probit_re', file = 'ProbitRE_ProbitRE.R', structure = 'panel',
  family = 'probit', has_selection = TRUE, parameters = c('sd_re_selection', 'sd_re_outcome', 'cor_re', 'cor_error'),
  likelihood = probit_re_probit_re_likelihood, panel = probit_re_probit_re_panel,
  outcome_rows = probit_re_rows, predictions = probit_re_probit_re_predictions,
  zero_correlations = list(cor_re = c('sd_re_selection', 'sd_re_outcome')))
