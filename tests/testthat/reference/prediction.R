# Posterior person effects by direct tensor integration, independent of the
# production posterior rules. Uses the base-R references in math.R.

# Posterior means and SDs of the person random intercept contributions
# sd_re_selection*u and sd_re_outcome*h given one panel's history. h = rho*u +
# sqrt(1 - rho^2)*v for selection models and h = u otherwise; u and v are
# independent standard normal coordinates (math.R ref_conditional_panel).
ref_posterior_effects = function(model, par, history, outer = 25L, inner = 35L) {
  spec = ref_model(model)
  rule = ref_normal_rule(outer)
  error_rule = ref_normal_rule(inner)
  grid = if (spec$selection) expand.grid(u = seq_len(outer), v = seq_len(outer)) else
    data.frame(u = seq_len(outer), v = rep(1L, outer))
  u = rule$nodes[grid$u]
  v = if (spec$selection) rule$nodes[grid$v] else rep(0, nrow(grid))
  weights = rule$weights[grid$u]
  if (spec$selection) weights = weights * rule$weights[grid$v]
  rho = ref_parameter(par, 'cor_re')
  h = if (spec$selection) rho*u + sqrt(1 - rho*rho)*v else u
  weights = weights * vapply(seq_len(nrow(grid)), function(k)
    ref_conditional_panel(model, par, history, u[k], v[k], error_rule)$value, numeric(1))
  weights = weights/sum(weights)
  moments = function(values) {
    mean = sum(weights*values)
    c(mean, sqrt(sum(weights*(values - mean)^2)))
  }
  out = stats::setNames(moments(ref_parameter(par, 'sd_re_outcome')*h), c('outcome_mean', 'outcome_sd'))
  if (spec$selection) out = c(stats::setNames(moments(ref_parameter(par, 'sd_re_selection')*u),
    c('selection_mean', 'selection_sd')), out)
  out
}
