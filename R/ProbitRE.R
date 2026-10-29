# Panel probit.
#
#   Y_it = 1{x_it'beta + o_it + lambda*v_i + eps_it > 0},  v_i, eps_it ~ N(0, 1) independent
#
#   symbol  parameter      meaning
#   beta    outcome.*      outcome coefficients
#   lambda  sd_re_outcome  person-level SD
#
# The model point m holds b_t = x_t'beta + o_t; at a standard-normal node z,
# v = z and B_t = b_t + lambda*v when the effect is active. Fixed GH is
# the following finite sum, evaluated in bounded batches:
#   for each node k: r_ik = sum_t log Phi(q_t B_tk),  q_t = 2y_t - 1
#   log L_i = log sum_k w_k exp(r_ik),  score_i = sum_k pi_ik g_i(z_k)  (pi: node posterior)
# Here g_i is the natural score of r_i at a fixed node. Read probit_re_rows
# for the observation partials and probit_re_panel for their parameter chains;
# model_assembly.R owns row sums and higher mixed derivatives.
#
# Routes (integrate_standalone, quadrature.R):
#   fixed GH                 integrate_index_block: probit_re_rows at every node, with
#                            node partials averaged before x is expanded
#   adaptive GH, history     integrate_person_effects; probit_re_rows supplies values
#                            and first scores through evaluate_standalone_nodes,
#                            while probit_re_panel supplies higher geometry derivatives
#   GL (moving cuts)         split_probit_panel (quadrature_split.R section 2)
# selection_component also uses this model for selection-only panel likelihoods,
# with the selection design and lambda = sd_re_selection. History retains nodes
# for prediction; its unreduced geometry is set by posterior_rule.

# log Phi(sign*B) with sign = 2y - 1, and its B derivatives.
probit_re_rows = function(B, y, log_factorial, m, order)
  local_univariate_probit(B, 2*y-1, 'b', order)

# r_i(z) = sum_t log Phi(q_t B_t), B_t = b_t + lambda*v, v = z;
# g_beta = sum_t x_t l_B, g_lambda = v sum_t l_B.
probit_re_panel = function(z, m, out) {
  e = person_effects(m, z, out, 'v')
  rows = node_rows(out, m$data$active)
  B = m$b[rows$row]+m$lambda*e$v[rows$unit]
  add_rows(out, e, probit_re_rows(B, rows$y, NULL, m, out$order), m$data$active,
    B = list(beta = 'x', lambda = e$v))
}

probit_re_predictions = function(context) {
  scale = if (context$conditional) 1 else sqrt(1+context$lambda^2)
  b = context$b/scale
  gradient_index = context$gradient_b/scale - if (context$conditional) context$zero else context$gradient_lambda*(b*context$lambda/scale^2)
  list(response = list(fit = stats::pnorm(b), gradient = gradient_index*stats::dnorm(b)),
       index = b, index_gradient = gradient_index, scale = scale)
}

# Natural layout and capabilities are declared beside the model equations.
probit_re_model = function() model_definition(
  id = 'probit_re', file = 'ProbitRE.R', structure = 'panel',
  family = 'probit', has_selection = FALSE, parameters = 'sd_re_outcome',
  likelihood = integrate_standalone, panel = probit_re_panel,
  outcome_rows = probit_re_rows, predictions = probit_re_predictions)
