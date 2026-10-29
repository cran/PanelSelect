# Panel Poisson.
#
#   Y_it | v_i ~ Poisson(exp(x_it'beta + o_it + lambda*v_i)),  v_i ~ N(0, 1)
#
#   symbol  parameter      meaning
#   beta    outcome.*      outcome coefficients
#   lambda  sd_re_outcome  person-level SD
#
# The model point m holds b_t = x_t'beta + o_t; at a standard-normal node z,
# v = z and B_t = b_t + lambda*v when the effect is active. Fixed GH is
# the following finite sum, evaluated in bounded batches:
#   for each node k: r_ik = sum_t [y_t B_tk - exp(B_tk) - log y_t!]
#   log L_i = log sum_k w_k exp(r_ik),  score_i = sum_k pi_ik g_i(z_k)  (pi: node posterior)
# Here g_i is the natural score of r_i at a fixed node. Read poisson_re_rows
# for the observation partials and poisson_re_panel for their parameter chains;
# model_assembly.R owns row sums and higher mixed derivatives.
#
# Routes (integrate_standalone, quadrature.R):
#   fixed GH                 integrate_index_block: poisson_re_rows at every node, with
#                            node partials averaged before x is expanded
#   adaptive GH, history     integrate_person_effects (quadrature.R)
#   GL (staged meshes)       fixed cells (quadrature_staged.R)
# In the latter routes, evaluate_standalone_nodes uses poisson_re_rows for
# values and first scores; poisson_re_panel supplies higher geometry derivatives.
# History retains nodes for prediction through posterior_rule.

poisson_re_rows = function(B, y, log_factorial, m, order) {
  mu = exp(B)
  # l(B) = y B - mu - log(y!); d l / d B = y - mu.
  out = local_block(y*B-mu-log_factorial, 'b', 1L, order, if (order) y-mu)
  if (order >= 2L) out$hessian[, 1L, 1L] = -mu
  if (order >= 3L) out$third[, 1L, 1L, 1L] = -mu
  out
}

# r_i(z) = sum_t log Pois(y_t; exp(B_t)), B_t = b_t + lambda*v, v = z;
# g_beta = sum_t x_t l_B, g_lambda = v sum_t l_B.
poisson_re_panel = function(z, m, out) {
  e = person_effects(m, z, out, 'v')
  rows = node_rows(out, m$data$active)
  B = m$b[rows$row]+m$lambda*e$v[rows$unit]
  add_rows(out, e, poisson_re_rows(B, rows$y, m$data$log_factorial[rows$row], m, out$order),
    m$data$active, B = list(beta = 'x', lambda = e$v))
}

count_prediction_mean = function(context, gamma = 0, gradient_gamma = context$zero) {
  logM = context$b+gamma^2/2+if (context$conditional) 0 else context$lambda^2/2
  gradient_log_M = context$gradient_b+gamma*gradient_gamma+
    if (context$conditional) context$zero else context$lambda*context$gradient_lambda
  M = exp(logM)
  list(response = list(fit = M, gradient = gradient_log_M*M), log_mean = logM,
       log_gradient = gradient_log_M)
}

poisson_re_predictions = function(context) count_prediction_mean(context)

# Natural layout and capabilities are declared beside the model equations.
poisson_re_model = function() model_definition(
  id = 'poisson_re', file = 'PoissonRE.R', structure = 'panel',
  family = 'poisson', has_selection = FALSE, parameters = 'sd_re_outcome',
  likelihood = integrate_standalone, panel = poisson_re_panel,
  outcome_rows = poisson_re_rows, predictions = poisson_re_predictions)
