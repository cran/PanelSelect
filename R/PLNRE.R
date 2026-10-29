# Panel Poisson-lognormal.
#
#   Y_it | v_i, e_it ~ Poisson(exp(x_it'beta + o_it + lambda*v_i + gamma*e_it)),
#   v_i, e_it ~ N(0, 1) independent: e_it is each row's own observation effect
#
#   symbol  parameter      meaning
#   beta    outcome.*      outcome coefficients
#   lambda  sd_re_outcome  person-level SD
#   gamma   sd_error       observation-level SD
#
# The model point m holds b_t = x_t'beta + o_t; at a standard-normal node z,
# v = z and B_t = b_t + lambda*v when the person effect is active. Each row
# integrates its own e before its contribution enters the person-level sum.
# With fixed GH at both levels, the calculation is:
#   for each node k: r_ik = sum_t log sum_j w_j Pois(y_t; exp(B_tk + gamma e_j))
#   log L_i = log sum_k w_k exp(r_ik),  score_i = sum_k pi_ik g_i(z_k)  (pi: node posterior)
# Here g_i is the natural score of r_i at fixed person nodes, after the row
# integral. Read pln_re_inner_kernel for its integrand and partials, pln_re_rows
# for the reduced row, and pln_re_panel for the person-effect parameter chains.
# quadrature_inner.R selects and reduces m$inner; model_assembly.R assembles
# person scores and higher mixed derivatives.
#
# Routes (integrate_standalone, quadrature.R):
#   fixed GH                 integrate_index_block: pln_re_rows at every node, with
#                            node partials averaged before x is expanded
#   adaptive GH, history     integrate_person_effects (quadrature.R)
#   GL (staged meshes)       fixed cells (quadrature_staged.R)
#   values and first scores  in those routes, pln_re_rows via evaluate_standalone_nodes;
#                            pln_re_panel supplies higher geometry derivatives
#   observation rule         fixed GH (default), adaptive GH (quadrature_inner_adaptive.R)
#                            or staged GL meshes (quadrature_staged_inner.R)
#   gamma fixed at 0         ordinary Poisson rows (poisson_re_rows)
# History retains person nodes for prediction through posterior_rule.

pln_re_inner_kernel = function(A, B, y, gamma, tau, e, log_weights, order, log_factorial,
                              parameters = TRUE) {
  eta = B+gamma*e; mu = exp(eta)
  coordinates = c('b', if (parameters) 'error')
  value = y*eta-mu-log_factorial+log_weights
  if (!order) return(local_block(value, coordinates, 1L, order))
  d_B = y-mu
  out = local_block(value, coordinates, 1L, order, if (parameters) cbind(d_B, e*d_B) else d_B)
  if (order >= 2L) {
    derivative = matrix(0, length(B), length(coordinates))
    derivative[, 1L] = 1
    if (parameters) derivative[, 2L] = e
    out = poisson_log_higher(out, derivative, mu, 1L, order)
  }
  out
}

pln_re_rows = function(B, y, log_factorial, m, order) {
  # Fixed gamma=0 removes the observation integral, leaving ordinary Poisson.
  if (is.null(m$inner)) return(poisson_re_rows(B, y, log_factorial, m, order))
  integrate_pln_inner(NULL, B, y, m$error, 0, order, log_factorial,
    m$inner, pln_re_inner_kernel, m$parameter_derivatives %||% TRUE)
}

# r_i(z) = sum_t log int Pois(y_t; exp(B_t + gamma e)) phi(e) de, B_t = b_t + lambda*v;
# g_beta = sum_t x_t l_B, g_lambda = v sum_t l_B, g_gamma = sum_t l_gamma.
pln_re_panel = function(z, m, out) {
  e = person_effects(m, z, out, 'v')
  rows = node_rows(out, m$data$active)
  B = m$b[rows$row]+m$lambda*e$v[rows$unit]
  add_rows(out, e, pln_re_rows(B, rows$y, m$data$log_factorial[rows$row], m, out$order),
    m$data$active, B = list(beta = 'x', lambda = e$v), error = list(gamma = 1))
}

pln_re_predictions = function(context)
  count_prediction_mean(context, context$error, context$gradient_error)

# Natural layout and capabilities are declared beside the model equations.
pln_re_model = function() model_definition(
  id = 'pln_re', file = 'PLNRE.R', structure = 'panel',
  family = 'pln', has_selection = FALSE, parameters = c('sd_re_outcome', 'sd_error'),
  likelihood = integrate_standalone, panel = pln_re_panel,
  outcome_rows = pln_re_rows, predictions = pln_re_predictions)
