# Panel probit selection with a Linear outcome.
#
#   D_it  = 1{w_it'alpha + o^s_it + delta*u_i + xi_it > 0}
#   Y*_it = x_it'beta + o^y_it + lambda*v_i + sigma*eps_it, observed when D_it = 1
#   (u_i, v_i) ~ N2(0, [1 rho; rho 1]),  (xi_it, eps_it) ~ N2(0, [1 tau; tau 1])
#
#   symbol  parameter        meaning
#   alpha   selection.*      selection coefficients
#   beta    outcome.*        outcome coefficients
#   delta   sd_re_selection  person-level selection SD
#   lambda  sd_re_outcome    person-level outcome SD
#   sigma   sd_error         observation-level outcome SD
#   rho     cor_re           correlation of the person effects
#   tau     cor_error        correlation of the observation shocks
#
# With exact reductions (the default), let y_obs contain this person's observed
# outcomes. The exact Gaussian marginal f_Y(y_obs) factors the likelihood as
#   L_i = f_Y(y_obs) * E[prod_t p_t(U,V) | y_obs].
# probit_re_linear_re_point supplies the conditioned person-effect maps.
# Given (U,V), A_t = a_t + delta*U, B_t = b_t + lambda*V and R_t = (y_t-B_t)/sigma,
#   p_t = Phi(-A_t) (d = 0),  Phi((A_t + tau R_t)/c_tau) (d = 1),
# where c_tau = sqrt(1-tau^2). Thus fixed GH of the conditional factor gives
#   r_ik = sum_t log p_t(U_ik,V_ik); log L_i = log f_Y(y_obs) + log sum_k w_k exp(r_ik).
# Without exact reductions, the selected rows also carry phi(R_t)/sigma and
# the person effects follow their prior maps. History uses this unreduced form
# to retain the latent variables needed by prediction (posterior_rule).
#
# Read linear_selection_log_contribution for observation partials, then
# probit_re_linear_re_panel for the natural first-score chains. The conditioning
# helpers below supply the exact marginal and maps; model_assembly.R adds the
# marginal and constructs latent/mixed derivatives from the same maps.
#
# Routes (probit_re_linear_re_likelihood):
#   fixed/adaptive GH            probit_re_linear_re_panel via integrate_person_effects
#   exact reductions, tau fixed 0, fixed GH
#                                probit_re_linear_re_fixed_block: selection moments
#                                averaged before the conditioned map is expanded
#   rho = tau = 0 fixed          ProbitRE selection factor plus the exact outcome factor
#                                (likelihood.R; adaptive GH keeps this model's geometry)
#   panels without outcomes      ProbitRE selection component (likelihood.R)
#   GL, tau fixed 0 (moving cuts) split_linear_panel (quadrature_split.R section 4)
#   other GL (staged meshes)     probit_re_linear_re_panel on fixed cells (quadrature_staged.R)
# Unselected rows use log Phi(-A), with d_A = -Mills(-A), without evaluating
# missing outcomes. Their local outcome partials are zero; under conditioning,
# the map of U can still carry outcome-parameter derivatives into this factor.

# Conditional row likelihood and natural first score -----------------------

linear_selection_log_contribution = function(A, B, y, sigma, tau, order, density) {
  n = length(A); coordinates = c('a', 'b', 'error', 'tau')
  R = (y-B)/sigma; c_tau = sqrt(1-tau*tau)
  zeta = (A+tau*R)/c_tau
  if (!order) {
    value = stats::pnorm(zeta, log.p = TRUE)
    if (density) value = value-R^2/2-log(2*pi)/2-log(sigma)
    return(local_block(value, coordinates, 2L, order))
  }
  # d log Phi(zeta) = Mills(zeta)*d zeta. The density below adds
  # R/sigma to d_B and (R^2-1)/sigma to d_sigma when it is retained.
  d_zeta_A = rep(1/c_tau, n)
  d_zeta_B = rep(-tau/(sigma*c_tau), n)
  d_zeta_sigma = -tau*R/(sigma*c_tau)
  d_zeta_tau = R/c_tau+zeta*tau/c_tau^2
  derivative = unname(cbind(d_zeta_A, d_zeta_B, d_zeta_sigma, d_zeta_tau))
  normal = normal_log_derivatives(zeta, order)
  # The gradient columns are d_A, d_B, d_sigma, d_tau of log Phi(zeta).
  out = local_block(normal$value, coordinates, 2L, order, derivative*normal$first)
  if (order >= 2L) {
    mixed = array(0, c(n, 2L, 4L))
    mixed[, 1L, 4L] = tau/c_tau^3
    mixed[, 2L, 3L] = tau/(sigma^2*c_tau)
    mixed[, 2L, 4L] = -1/(sigma*c_tau^3)
    out = normal_log_higher(out, normal, derivative, mixed, 2L, order)
  }
  if (!density) return(out)
  out$value = out$value - R^2/2 - log(2*pi)/2 - log(sigma)
  out$gradient[, 2L] = out$gradient[, 2L] + R/sigma
  out$gradient[, 3L] = out$gradient[, 3L] + (R^2-1)/sigma
  if (order >= 2L) {
    out$hessian[, 2L, 2L] = out$hessian[, 2L, 2L] - 1/sigma^2
    out$hessian[, 2L, 3L] = out$hessian[, 2L, 3L] - 2*R/sigma^2
  }
  if (order >= 3L) out$third[, 2L, 2L, 3L] = out$third[, 2L, 2L, 3L] + 2/sigma^3
  out
}

# Panel kernel. With exact reductions (the default) the outcome density is the
# exact marginal in m$marginal and u, v are the outcome-conditioned maps of
# linear_conditional_maps(), which depend on beta,
# lambda, sigma and rho. Otherwise u = z1, v = rho*z1 + sqrt(1-rho^2)*z2 and the
# selected rows keep the outcome density. The shared evaluator includes the
# marginal at each node; integration may pull it outside the node sum once.
probit_re_linear_re_panel = function(z, m, out) {
  if (m$gaussian && m$tau_zero) {
    # Exact outcome and fixed tau = 0: r_i = sum_t log Phi(s_t A_t), s_t = 2d_t - 1,
    # with U conditioned on the observed outcomes.
    e = person_effects(m, z, out, 'u')
    rows = node_rows(out, m$data$active)
    A = m$a[rows$row]+m$delta*e$u[rows$unit]
    return(add_rows(out, e, local_univariate_probit(A, 2*m$d[rows$row]-1, 'a', out$order),
      m$data$active, A = list(alpha = 'w', u = m$delta*e$du, delta = e$u)))
  }
  e = person_effects(m, z, out)
  # dA = w'dalpha + delta dU + U ddelta and dB = x'dbeta + lambda dV + V dlambda,
  # where dU and dV carry the conditioning columns. Coincident lambda terms of B
  # are combined before multiplying by sum_t l_B.
  dA = list(alpha = 'w', u = m$delta*e$du, delta = e$u)
  dB_v = e$dv*m$lambda
  at = match(out$columns$lambda, e$columns$v)
  if (length(at) && !is.na(at)) dB_v[, at] = dB_v[, at]+e$v
  dB = list(beta = 'x', v = dB_v, lambda = if (!length(at) || is.na(at)) e$v)
  unselected = node_rows(out, m$data$unselected)
  A = m$a[unselected$row]+m$delta*e$u[unselected$unit]
  out = add_rows(out, e, local_univariate_probit(A, -1, 'a', out$order),
    m$data$unselected, A = dA)
  selected = node_rows(out, m$data$selected)
  A = m$a[selected$row]+m$delta*e$u[selected$unit]
  B = m$b[selected$row]+m$lambda*e$v[selected$unit]
  # log Phi((A + tau R)/c_tau), R = (y - B)/sigma, plus log phi(R) - log sigma
  # unless the outcome density is already in m$marginal.
  add_rows(out, e, linear_selection_log_contribution(A, B, selected$y, m$error, m$tau,
      out$order, density = !m$gaussian), m$data$selected,
    A = dA, B = dB, error = list(sigma = 1), tau = list(tau = 1))
}

# Likelihood: all Linear routes. With exact reductions and tau fixed at 0, the
# selection factor given the outcomes is a probit in U alone: with rho also fixed
# at 0, U is independent of y and the factor is the standalone probit of the
# selection component (ProbitRE.R; adaptive GH keeps this model's geometry); with
# fixed GH and one effect, probit_re_linear_re_fixed_block averages selection
# moments before expanding the conditioned map. The exact outcome marginal is
# added once per unit.
probit_re_linear_re_likelihood = function(theta, prepared, control, fixed, q, score, posterior) {
  m = make_model_context(theta, prepared, control, fixed)
  adaptive = q > 0 && control$quadrature == 'agh'
  if (!posterior && m$gaussian && m$tau_zero) {
    if (is_fixed_zero('cor_re', theta, fixed)) {
      selection = prepared$integration_cache$selection %||%
        prepare_model_data(selection_component(prepared))
      return(integrate_person_effects(theta, prepared, control, fixed, q, score, m = m,
        block = if (!adaptive) integrate_index_block,
        kernel_point = make_model_context(theta, selection, control, fixed),
        kernel_data = selection))
    }
    if (!adaptive && q == 1L)
      return(integrate_person_effects(theta, prepared, control, fixed, q, score, m = m,
        block = probit_re_linear_re_fixed_block))
  }
  integrate_person_effects(theta, prepared, control, fixed, q, score, posterior, m = m)
}

# With rho and tau fixed at 0 the outcomes are independent of selection; their
# exact factor is the marginal outcome likelihood (likelihood.R).
probit_re_linear_re_outcome_factor = function(theta, prepared)
  linear_outcome_marginal(theta, prepared, conditional = FALSE)

# Exact outcome conditioning. With exact reductions (the default) phi(y | b,
# lambda, sigma) is integrated analytically: m$marginal holds its log value and
# natural gradient per unit, and u, v become the outcome-conditioned maps
# V = mu_v + s_v z_v, U = rho V + c_rho z_u (tau free) or
# U = rho mu_v + sqrt(1 - rho^2 + rho^2 s_v^2) z (tau fixed 0), whose locations
# and loadings depend on beta, lambda, sigma and rho.
probit_re_linear_re_point = function(m, theta, prepared) {
  if (!m$gaussian) return(standard_person_effects(m, theta, prepared))
  outcome = linear_outcome_marginal(theta, prepared, predictor = m$b)
  m$marginal = outcome[c('value', 'gradient')]
  maps = linear_conditional_maps(outcome, m$u, m$v, m$rho, m$tau_zero,
    m$has_u, m$has_v, m$rho_index)
  m$u = maps$u; m$v = maps$v
  m
}

# Fixed tau=0: U = mu_u + sqrt(k_u)*z; average selection moments before
# expanding the derivatives of mu_u and sqrt(k_u) across natural parameters.
# The shared integrator adds the exact outcome marginal once after all blocks.
probit_re_linear_re_fixed_block = function(theta, prepared, context, base, log_weights, score) {
  block = context$data$active; count = length(block$rows)
  n = prepared$n_panels; p = length(theta); size = nrow(base)
  z = base[, 1L]
  U = outer(context$u$loading[, 1L], z)+context$u$location
  A = matrix(context$a[block$rows], count, size)+
    context$delta*U[block$group, , drop = FALSE]
  sign = 2*prepared$d[block$rows]-1
  local = normal_log_derivatives(as.vector(A)*sign, as.integer(score))
  value = group_sum(matrix(local$value, count, size), block$ptr)
  mixture = logsum_weights(as.vector(value)+rep(log_weights, each = n), n)
  gradient = matrix(0, n, p)
  if (score) {
    d_A = matrix(local$first, count, size)*sign
    total = group_sum(d_A, block$ptr)
    total[mixture$weights == 0] = 0
    weighted = total*mixture$weights
    # S = sum_t d_A: dU contributes delta*(dmu_u E[S] + dsqrt(k_u) E[z*S]).
    gradient = context$delta*(context$u$location_gradient*rowSums(weighted)+
      matrix(context$u$loading_gradient[, 1L, ], n, p)*drop(weighted %*% z))
    columns = context$score_columns
    gradient[, columns$delta] = gradient[, columns$delta]+rowSums(weighted*U)
    products = group_weighted_design(block$w, d_A, block$ptr, mixture$weights)
    if (length(columns$alpha)) gradient[, columns$alpha] =
      gradient[, columns$alpha, drop = FALSE]+products[, seq_along(columns$alpha), drop = FALSE]
  }
  list(loglik = mixture$value, score = gradient)
}

# Exact outcome marginal and differentiated conditioning -------------------

linear_outcome_marginal = function(theta, prepared, conditional = TRUE, predictor = NULL) {
  prepared = prepare_model_data(prepared)
  data = prepared$model_data; observed = data$selected
  n = prepared$n_panels; p = length(theta); m = observed$count
  lambda_name = prepared$outcome_scale %||% 'sd_re_outcome'
  relevant = c(data$outcome_names, lambda_name, 'sd_error')
  key = list(parameters = theta[relevant], layout = names(theta), conditional = conditional)
  memo = prepared$integration_cache$memo
  if (!is.null(memo) && identical(memo$gaussian_key, key)) return(memo$gaussian_value)

  ix = match(data$outcome_names, names(theta))
  il = match(lambda_name, names(theta)); ie = match('sd_error', names(theta))
  lambda = unname(theta[il]); sigma = unname(theta[ie])
  if (is.null(predictor)) predictor = (if (length(ix))
    drop(prepared$x %*% theta[ix]) else rep(0, prepared$n_rows))+prepared$offset_outcome
  residual = observed$y-predictor[observed$rows]
  total = group_sum(cbind(residual), observed$ptr)[, 1L]
  # Within-panel centring: with W the centred sum of squares, the quadratic form
  # is W/sigma^2 + T^2/(mV). The raw form cancels a sum of squares against the
  # squared total when m*lambda^2 dominates sigma^2.
  count = pmax(m, 1)
  rows = rep(seq_len(n), m)
  centered = residual-total[rows]/count[rows]
  # Centring x as well keeps the panel-mean roundoff of the residuals out of
  # the coefficient score; an intercept column contributes exactly zero here.
  x_centered = observed$x-observed$sum_x[rows, , drop = FALSE]/count[rows]
  sums = group_sum(cbind(centered^2, x_centered*centered), observed$ptr)
  within = sums[, 1L]
  variance = sigma^2+m*lambda^2
  value = -.5*(m*log(2*pi)+(m-1)*log(sigma^2)+log(variance)+
                 within/sigma^2+total^2/(count*variance))
  gradient = matrix(0, n, p)
  if (length(ix)) gradient[, ix] = sums[, -1L, drop = FALSE]/sigma^2 +
    observed$sum_x*(total/(count*variance))
  gradient[, il] = lambda*(total^2/variance^2-m/variance)
  gradient[, ie] = -(m-1)/sigma-sigma/variance+within/sigma^3+
    sigma*total^2/(count*variance^2)
  value[m == 0] = 0; gradient[m == 0, ] = 0
  result = list(value = value, gradient = gradient)
  if (conditional) {
    mean_v = lambda*total/variance; var_v = sigma^2/variance; scale_v = sqrt(var_v)
    dmean = dvar = matrix(0, n, p)
    if (length(ix)) dmean[, ix] = -observed$sum_x*(lambda/variance)
    dmean[, il] = total*(sigma^2-m*lambda^2)/variance^2
    dmean[, ie] = -2*sigma*lambda*total/variance^2
    dvar[, il] = -2*m*lambda*sigma^2/variance^2
    dvar[, ie] = 2*m*sigma*lambda^2/variance^2
    result = c(result, list(mean = mean_v, variance = var_v, scale = scale_v,
      mean_gradient = dmean, variance_gradient = dvar, scale_gradient = dvar/(2*scale_v)))
  }
  if (!is.null(memo)) { memo$gaussian_key = key; memo$gaussian_value = result }
  result
}

# With tau fixed at 0, U given the observed outcomes is normal with mean
# rho*mu_v and standard deviation sqrt(1 - rho^2 + rho^2 s_v^2) (V integrated out
# of the selection factor). Location and scale with their natural-parameter
# derivatives; rows are units (one row per panel for moving-cut GL).
linear_selection_given_outcome = function(rho, mean_v, var_v, dmean, dvar, rho_index) {
  scale = sqrt(1-rho^2+rho^2*var_v)
  location_gradient = rho*dmean
  scale_gradient = rho^2*dvar/(2*scale)
  if (!is.na(rho_index)) {
    location_gradient[, rho_index] = location_gradient[, rho_index]+mean_v
    scale_gradient[, rho_index] = scale_gradient[, rho_index]+rho*(var_v-1)/scale
  }
  list(location = rho*mean_v, scale = scale, location_gradient = location_gradient,
       scale_gradient = scale_gradient)
}

# Exact Y-observed conditioning: V = mean_v + scale_v*z2 and
# U = rho*V + sqrt(1-rho^2)*z1. Fixed tau=0 integrates V out of selection.
# The supplied empty maps retain the shared layout and all natural columns.
linear_conditional_maps = function(outcome, u, v, rho, tau_zero, has_u, has_v, rho_index) {
  n = length(outcome$value); q = ncol(u$loading); ir = rho_index
  mean_v = outcome$mean; var_v = outcome$variance; scale_v = outcome$scale
  dmean = outcome$mean_gradient; dvar = outcome$variance_gradient; dscale = outcome$scale_gradient
  if (tau_zero) {
    v$location = mean_v; v$location_gradient = dmean
    if (has_u) {
      given = linear_selection_given_outcome(rho, mean_v, var_v, dmean, dvar, ir)
      u$location = given$location; u$location_gradient = given$location_gradient
      u$loading[, 1L] = given$scale
      u$loading_gradient[, 1L, ] = given$scale_gradient
    }
  } else {
    if (has_v) {
      k = if (has_u) 2L else 1L
      v$location = mean_v; v$location_gradient = dmean
      v$loading[, k] = scale_v; v$loading_gradient[, k, ] = dscale
    }
    if (has_u) {
      u$location = rho*v$location; u$location_gradient = rho*v$location_gradient
      u$loading = rho*v$loading; u$loading[, 1L] = sqrt(1-rho^2)
      u$loading_gradient = rho*v$loading_gradient
      if (!is.na(ir)) {
        u$location_gradient[, ir] = u$location_gradient[, ir]+v$location
        u$loading_gradient[, , ir] = matrix(u$loading_gradient[, , ir], n, q)+v$loading
        u$loading_gradient[, 1L, ir] = -rho/sqrt(1-rho^2)
      }
    }
  }
  list(u = u, v = v)
}

probit_re_linear_re_predictions = function(context) {
  M = context$b; gradient_M = context$gradient_b
  response = list(fit = M, gradient = gradient_M)
  if (!context$joint_needed) return(list(response = response))
  covariance = if (context$conditional) context$error*context$tau else
    (context$delta*context$lambda*context$rho+context$error*context$tau)/context$selection_scale
  gradient_covariance = if (context$conditional) context$tau*context$gradient_error+
      context$error*context$gradient_tau else
    (context$lambda*context$rho*context$gradient_delta+
      context$delta*context$rho*context$gradient_lambda+
      context$delta*context$lambda*context$gradient_rho+
      context$tau*context$gradient_error+
      context$error*context$gradient_tau)/
      context$selection_scale-context$gradient_delta*(covariance*context$delta/context$selection_scale^2)
  mills = prediction_mills(context$a); density = stats::dnorm(context$a)
  C = M+covariance*mills
  gradient_C = gradient_M+gradient_covariance*mills-context$gradient_a*(covariance*mills*(context$a+mills))
  J = M*context$selection$fit+covariance*density
  gradient_J = gradient_M*context$selection$fit+
      context$selection$gradient*M+gradient_covariance*density-
      context$gradient_a*(covariance*context$a*density)
  list(response = response, joint = list(fit = J, gradient = gradient_J),
       selected = list(fit = C, gradient = gradient_C))
}

# Natural layout and capabilities are declared beside the model equations.
probit_re_linear_re_model = function() model_definition(
  id = 'probit_re_linear_re', file = 'ProbitRE_LinearRE.R', structure = 'panel',
  family = 'gaussian', has_selection = TRUE, parameters = c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error'),
  likelihood = probit_re_linear_re_likelihood, panel = probit_re_linear_re_panel,
  point = probit_re_linear_re_point, exact_outcome = probit_re_linear_re_outcome_factor,
  predictions = probit_re_linear_re_predictions,
  positive_scales = 'sd_error',
  zero_correlations = list(cor_re = c('sd_re_selection', 'sd_re_outcome')))
