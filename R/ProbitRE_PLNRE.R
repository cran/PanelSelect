# Panel probit selection with a Poisson-lognormal outcome.
#
#   D_it = 1{w_it'alpha + o^s_it + delta*u_i + xi_it > 0}
#   Y*_it | v_i, e_it ~ Poisson(exp(x_it'beta + o^y_it + lambda*v_i + gamma*e_it)),
#   observed when D_it = 1; (u_i, v_i) ~ N2(0, [1 rho; rho 1]),
#   (xi_it, e_it) ~ N2(0, [1 tau; tau 1]): each row has its own e, shared only
#   within its selection and count contribution
#
#   symbol  parameter        meaning
#   alpha   selection.*      selection coefficients
#   beta    outcome.*        outcome coefficients
#   delta   sd_re_selection  person-level selection SD
#   lambda  sd_re_outcome    person-level outcome SD
#   gamma   sd_error         observation-level SD
#   rho     cor_re           correlation of the person effects
#   tau     cor_error        correlation of xi and e
#
# The model point m holds a_t and b_t. With both person effects active, normal
# nodes z give u = z1, v = rho*z1 + sqrt(1 - rho^2)*z2, so A_t = a_t + delta*u and
# B_t = b_t + lambda*v. Each selected row integrates its e by the observation
# rule m$inner (quadrature_inner.R). With fixed GH at both levels:
#   for each node k: r_ik = sum_{d=0} log Phi(-A_tk) + sum_{d=1} log p_t(A_tk, B_tk),
#     p_t = sum_j w_j Phi((A + tau e_j)/c_tau) Pois(y_t; exp(B + gamma e_j))
#   log L_i = log sum_k w_k exp(r_ik),  score_i = sum_k pi_ik g_i(z_k)  (pi: node posterior)
# c_tau = sqrt(1-tau^2); g_i is the natural score at fixed person nodes after
# each row integral. Read probit_re_pln_re_inner_kernel for the integrand and
# partials, probit_re_pln_re_rows for its integral, and probit_re_pln_re_panel
# for the observed likelihood and parameter chains. Unselected rows use only
# log Phi(-A). Shared integration reduces observation nodes before person nodes.
#
# Routes (likelihood.R dispatches reductions before this model's entry):
#   fixed/adaptive GH, history   probit_re_pln_re_panel via integrate_person_effects
#   fixed GH, tau fixed 0        integrate_selection_fixed: Phi(s_t A_t) reused across v | u
#   rho = tau = 0 fixed          independent ProbitRE and PLNRE factors (likelihood.R)
#   panels without outcomes      ProbitRE selection component (likelihood.R)
#   GL (staged meshes)           probit_re_pln_re_panel on fixed cells (quadrature_staged.R)
#   GL at both levels, tau fixed 0
#                                physical (u,v) rule when eligible (quadrature_staged_pln.R)
#   observation rule             fixed GH (default), adaptive GH (quadrature_inner_adaptive.R)
#                                or staged GL meshes (quadrature_staged_inner.R)
# GH factorization and nested shortcuts require exact_reductions; the nested
# fixed rule also requires both person effects active and no history output.
# A free tau at zero still needs its coupled derivative. History retains the
# unreduced person-effect representation and does not support GL rules.

probit_re_pln_re_inner_kernel = function(A, B, y, gamma, tau, e, log_weights, order, log_factorial,
    parameters = TRUE, normal_evaluate = normal_log_derivatives) {
  eta = B+gamma*e; mu = exp(eta)
  c_tau = sqrt(1-tau*tau); zeta = (A+tau*e)/c_tau
  coordinates = c('a', 'b', if (parameters) c('error', 'tau')); n = length(B)
  # Preserve Poisson + log GH weight + log selection in that order.
  log_poisson = y*eta-mu-log_factorial+log_weights
  if (!order) return(local_block(log_poisson+stats::pnorm(zeta, log.p = TRUE),
    coordinates, 2L, order))
  normal = normal_evaluate(zeta, order)
  d_A = (1/c_tau)*normal$first
  d_B = y-mu
  if (parameters) {
    d_gamma = e*d_B
    d_zeta_tau = e/c_tau+zeta*tau/c_tau^2
    d_tau = d_zeta_tau*normal$first
  }
  # Separate compact components preserve zero*infinite trial behavior and the
  # original first/higher addition order; impossible inner mass is masked later.
  count = local_block(log_poisson, coordinates, 2L, order,
    if (parameters) cbind(0*d_B, d_B, d_gamma, 0*d_B) else cbind(0*d_B, d_B))
  selection = local_block(normal$value, coordinates, 2L, order,
    if (parameters) cbind(d_A, 0*normal$first, 0*normal$first, d_tau) else cbind(d_A, 0*normal$first))
  if (order >= 2L) {
    derivative = matrix(0, n, length(coordinates))
    derivative[, 2L] = 1
    if (parameters) derivative[, 3L] = e
    count = poisson_log_higher(count, derivative, mu, 2L, order)
    derivative[] = 0; derivative[, 1L] = 1/c_tau
    if (parameters) derivative[, 4L] = d_zeta_tau
    mixed = array(0, c(n, 2L, length(coordinates)))
    if (parameters) mixed[, 1L, 4L] = tau/c_tau^3
    selection = normal_log_higher(selection, normal, derivative, mixed, 2L, order)
  }
  local_add(count, selection)
}

# Selected row with free tau: log int Phi((A + tau e)/c_tau) Pois(y; e^{B + gamma e})
# phi(e) de by the observation rule m$inner (fixed GH, adaptive GH or staged GL).
probit_re_pln_re_rows = function(A, B, y, log_factorial, m, order) {
  kernel = probit_re_pln_re_inner_kernel
  if (!is.null(m$selection_normal)) kernel = function(...)
    probit_re_pln_re_inner_kernel(..., normal_evaluate = m$selection_normal)
  integrate_pln_inner(A, B, y, m$error, m$tau, order, log_factorial,
    m$inner, kernel, m$parameter_derivatives %||% TRUE)
}

# r_i(z) = sum_{d=0} log Phi(-A_t) + sum_{d=1} log p_t(A_t, B_t), A_t = a_t + delta*u,
# B_t = b_t + lambda*v. With tau or gamma fixed at 0, log p_t = log Phi(A_t) +
# log PLN(y_t | B_t) (or log Pois when gamma = 0).
probit_re_pln_re_panel = function(z, m, out) {
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
  log_factorial = m$data$log_factorial[selected$row]
  if (!m$tau_zero && !is.null(m$inner))
    return(add_rows(out, e, probit_re_pln_re_rows(A, B, selected$y, log_factorial, m, out$order),
      m$data$selected, A = dA, B = dB, error = list(gamma = 1), tau = list(tau = 1)))
  out = add_rows(out, e, local_univariate_probit(A, 1, 'a', out$order), m$data$selected, A = dA)
  add_rows(out, e, pln_re_rows(B, selected$y, log_factorial, m, out$order), m$data$selected,
    B = dB, error = list(gamma = 1))
}

# Likelihood. With tau fixed at 0 and fixed GH over both effects, Phi(s_t A_t)
# depends on u only, so the rule conditions on u and integrates v | u over the
# outcome nodes (integrate_selection_fixed). With free tau, the fixed tensor
# revisits identical A_t = a_t + delta*z1 across outcome nodes, so the selection
# normal derivatives are cached within this evaluation. Otherwise
# integrate_person_effects integrates the panel kernel.
probit_re_pln_re_likelihood = function(theta, prepared, control, fixed, q, score, posterior) {
  m = make_model_context(theta, prepared, control, fixed)
  fixed_tensor = q == 2L && control$quadrature != 'agh' && !posterior
  if (fixed_tensor && isTRUE(control$exact_reductions) && m$tau_zero)
    return(integrate_selection_fixed(theta, prepared, m, control$.rules$outer[[3L]], score))
  if (fixed_tensor && score && !m$tau_zero && !is.null(m$inner))
    m$selection_normal = cache_normal_derivatives()
  integrate_person_effects(theta, prepared, control, fixed, q, score, posterior, m = m)
}

probit_re_pln_re_predictions = function(context) {
  gamma = context$error; gradient_gamma = context$gradient_error
  mean = count_prediction_mean(context, gamma, gradient_gamma)
  if (!context$joint_needed) return(mean)
  z = context$a+if (context$conditional) context$tau*gamma else
    (context$delta*context$lambda*context$rho+context$tau*gamma)/context$selection_scale
  gradient_z = context$gradient_a+if (context$conditional) context$tau*gradient_gamma+
      gamma*context$gradient_tau else
    (context$lambda*context$rho*context$gradient_delta+
      context$delta*context$rho*context$gradient_lambda+
      context$delta*context$lambda*context$gradient_rho+
      context$tau*gradient_gamma+
      gamma*context$gradient_tau)/
      context$selection_scale-
      context$gradient_delta*((context$delta*context$lambda*context$rho+
      context$tau*gamma)*context$delta/context$selection_scale^3)
  logJ = mean$log_mean+stats::pnorm(z, log.p = TRUE)
  gradient_log_J = mean$log_gradient+gradient_z*prediction_mills(z)
  J = exp(logJ); C = exp(logJ-context$selection$log_fit)
  list(response = mean$response, joint = list(fit = J, gradient = gradient_log_J*J),
    selected = list(fit = C, gradient = (gradient_log_J-context$selection$log_gradient)*C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_re_pln_re_model = function() model_definition(
  id = 'probit_re_pln_re', file = 'ProbitRE_PLNRE.R', structure = 'panel',
  family = 'pln', has_selection = TRUE, parameters = c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error'),
  likelihood = probit_re_pln_re_likelihood, panel = probit_re_pln_re_panel,
  outcome_rows = pln_re_rows, predictions = probit_re_pln_re_predictions,
  zero_correlations = list(cor_re = c('sd_re_selection', 'sd_re_outcome'), cor_error = 'sd_error'))
