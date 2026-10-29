# One likelihood evaluation: the log likelihood and natural score of every
# independent unit at theta.
#
#   route                                  calculation
#   staged GL meshes                       evaluate_staged_integration (quadrature_staged.R)
#   GL with moving cuts                    evaluate_split_binary (quadrature_split.R)
#   cross-sectional rows                   the model's likelihood entry, no person effects
#   rho, tau fixed at 0 with fixed GH      integrate_independent: selection x outcome factors
#   panels without observed outcomes       their selection factor (ProbitRE), then the rest
#   rho, tau fixed at 0 with adaptive GH   integrate_independent, unless the model integrates
#                                          an exact outcome factor itself (Linear)
#   otherwise                              integrate_prepared: the model's likelihood entry
#
# This is the authoritative route order; model files describe their local
# calculations and optional reductions. Reductions require exact_reductions
# where checked below. A correlation fixed at zero can change a route; a free
# correlation evaluated at zero retains its coupled parameter derivative.
# Returned scores use natural parameter columns with fixed columns masked to
# zero. The GH routes mask them here; the GL evaluators mask their own results.
# Lower-level kernels can still calculate derivatives of fixed parameters.

evaluate_likelihood = function(theta, prepared, control, fixed = numeric(), score = TRUE) {
  control = current_control(control)
  stop_removed_settings(control)
  if (uses_staged_integration(prepared, control, fixed))
    return(evaluate_staged_integration(theta, prepared, control, fixed, score))
  if (uses_gl_quadrature(control))
    return(evaluate_split_binary(theta, prepared, control, fixed, score))
  prepared = prepare_integration(prepared, control)
  control = prepared$integration_cache$control
  # With rho fixed at 0 and tau fixed at 0 (or absent, as for Poisson outcomes),
  # selection and outcome are independent factors of each panel's likelihood.
  independent = prepared$has_selection && isTRUE(control$exact_reductions) &&
    is_fixed_zero('cor_re', theta, fixed) && is_fixed_zero('cor_error', theta, fixed)
  # The remaining routes of the table, in order. Cross-sectional observations
  # have no panels to split or factor; a panel without observed outcomes is
  # exactly a standalone selection panel.
  result = if (identical(prepared$structure, 'cross_sectional')) {
    integrate_prepared(theta, prepared, control, fixed, score)
  } else if (independent && control$quadrature == 'gh') {
    integrate_independent(theta, prepared, control, fixed, score)
  } else if (prepared$has_selection && any(!prepared$integration_cache$observed)) {
    integrate_unobserved(theta, prepared, control, fixed, score)
  } else if (independent && is.null(prepared_model(prepared)$exact_outcome)) {
    integrate_independent(theta, prepared, control, fixed, score)
  } else integrate_prepared(theta, prepared, control, fixed, score)
  if (length(fixed)) result$score[, intersect(names(fixed), names(theta))] = 0
  result
}

# Compatibility name used by existing callers.
evaluate_model = evaluate_likelihood

prepared_model = function(prepared)
  resolve_model(prepared$family, prepared$has_selection, prepared$structure, prepared$model_id)

# Panels without an observed outcome contribute their selection factor (the
# ProbitRE selection component); the other panels follow every route again.
integrate_unobserved = function(theta, prepared, control, fixed, score) {
  cache = prepared$integration_cache; observed = cache$observed
  result = list(loglik = numeric(prepared$n_panels),
    score = matrix(0, prepared$n_panels, length(theta), dimnames = list(NULL, names(theta))))
  no = which(!observed); yes = which(observed)
  a = integrate_component(theta, cache$selection_only, control, fixed, score)
  result$loglik[no] = a$loglik; result$score[no, ] = a$score
  if (length(yes)) {
    b = evaluate_likelihood(theta, cache$with_outcomes, control, fixed, score)
    result$loglik[yes] = b$loglik; result$score[yes, ] = b$score
  }
  result$diagnostics = list(method = control$quadrature,
    selection_only_panels = prepared$id_values[no],
    observed_panels = prepared$id_values[yes], selection_only = a$diagnostics,
    observed = if (length(yes)) b$diagnostics else NULL, order_checked = FALSE)
  result
}

# Independent factors: the selection component (ProbitRE) on every panel, plus
# the outcome factor on panels with observed outcomes. The outcome factor is the
# model's exact outcome likelihood when it has one (Linear), otherwise its
# standalone outcome model (ProbitRE, PoissonRE or PLNRE).
integrate_independent = function(theta, prepared, control, fixed, score) {
  cache = prepared$integration_cache; observed = cache$observed
  result = integrate_component(theta, cache$selection, control, fixed, score)
  selection_diagnostics = result$diagnostics
  outcome_diagnostics = NULL
  if (any(observed)) {
    exact_outcome = prepared_model(prepared)$exact_outcome
    if (!is.null(exact_outcome)) {
      marginal = exact_outcome(theta, prepared)
      result$loglik = result$loglik + marginal$value
      if (score) result$score = result$score + marginal$gradient
      outcome_diagnostics = selection_diagnostics
    } else {
      part = if (all(observed)) prepared else cache$with_outcomes
      outcome = part$integration_cache$outcome
      value = integrate_component(theta, outcome, outcome$integration_cache$control, fixed, score)
      result$loglik[observed] = result$loglik[observed] + value$loglik
      result$score[observed, ] = result$score[observed, , drop = FALSE] + value$score
      outcome_diagnostics = list(method = control$quadrature, factorized = TRUE,
        selection = selection_diagnostics, outcome = value$diagnostics, order_checked = FALSE)
    }
  }
  result$diagnostics = if (all(observed)) outcome_diagnostics else
    list(method = control$quadrature, selection_only_panels = prepared$id_values[!observed],
      observed_panels = prepared$id_values[observed], selection_only = selection_diagnostics,
      observed = outcome_diagnostics, order_checked = FALSE)
  if (any(!is.finite(result$loglik)) || (score && any(!is.finite(result$score))))
    integration_error('Likelihood or score is not finite at this parameter vector')
  result
}

# One component (a selection or outcome factor) on the parameters it depends
# on, its score scattered back to every natural parameter. The PLN outcome
# factor also depends on sd_error.
integrate_component = function(theta, prepared, control, fixed, score) {
  relevant = c(prepared$model_data$outcome_names,
            prepared$outcome_scale %||% 'sd_re_outcome',
            if (prepared$family == 'pln') 'sd_error')
  take = names(theta) %in% relevant
  # Independent equations often keep the same coordinates across information
  # calculations. This private per-fit cache never substitutes a different point.
  key = list(theta = theta[take], fixed = fixed[names(fixed) %in% names(theta)[take]],
             score = score)
  memo = prepared$integration_cache$memo
  if (!is.null(memo) && identical(memo$key, key)) component = memo$value else {
    component = integrate_prepared(key$theta, prepared, control, key$fixed, score)
    if (!is.null(memo)) { memo$key = key; memo$value = component }
  }
  result = matrix(0, prepared$n_panels, length(theta), dimnames = list(NULL, names(theta)))
  result[, take] = component$score
  component$score = result
  component
}

# Integration caches ---------------------------------------------------------
# Built once per evaluator and settings: model data, rules, a memo and the
# selection, selection-only, observed-outcome and outcome components.

# Staged recorders, clocks and inner budgets are fresh on every evaluation but
# leave the rules unchanged; a cache hit re-attaches the current ones.
staged_transients = c('.staged_recorder', '.staged_clock', '.staged_inner_budget')

attach_staged_transients = function(cache, control) {
  for (field in staged_transients) cache$control[[field]] = control[[field]]
  for (part in c('selection', 'selection_only', 'with_outcomes', 'outcome'))
    if (!is.null(cache[[part]])) cache[[part]]$integration_cache =
      attach_staged_transients(cache[[part]]$integration_cache, control)
  cache
}

prepare_integration = function(prepared, control) {
  settings = control
  settings[c('.rules', staged_transients)] = NULL
  cached = prepared$integration_cache
  if (!is.null(cached) && identical(cached$settings, settings)) {
    if (any(staged_transients %in% c(names(control), names(cached$control))))
      prepared$integration_cache = attach_staged_transients(cached, control)
    return(prepared)
  }
  prepared = prepare_model_data(prepared)
  control = prepare_rule_control(control,
    if (identical(prepared$structure, 'cross_sectional')) 0L else 2L)
  cache = list(settings = settings, control = control, memo = new.env(parent = emptyenv()))
  if (prepared$has_selection && !identical(prepared$structure, 'cross_sectional')) {
    observed = tabulate(prepared$group[prepared$d == 1], prepared$n_panels) > 0
    cache$observed = observed
    cache$selection = prepare_integration(selection_component(prepared), control)
    if (any(!observed)) {
      cache$selection_only = prepare_integration(
        selection_component(subset_panels(prepared, which(!observed))), control)
      if (any(observed)) cache$with_outcomes = prepare_integration(
        subset_panels(prepared, which(observed)), control)
    }
    # The standalone outcome model of the independent factorization; a model
    # with an exact outcome likelihood (Linear) needs none.
    if (all(observed) && is.null(prepared_model(prepared)$exact_outcome)) {
      outcome = prepared
      outcome$has_selection = FALSE
      outcome$model_id = NULL
      outcome$active = prepared$d == 1
      outcome$model_data = NULL
      outcome$integration_cache = NULL
      c_outcome = control
      if (length(control$nodes) > 1L) c_outcome$nodes = control$nodes[2L]
      cache$outcome = prepare_integration(outcome, c_outcome)
    }
  }
  prepared$integration_cache = cache
  prepared
}

# Attach one private transaction workspace to the already prepared components.
# It belongs to an evaluator, never to the fitted object's prepared data. Paths
# distinguish joint, selection-only and independent factors without hashing IDs
# or copying mode matrices on every evaluation.
prepare_mode_reuse = function(prepared) {
  state = new.env(parent = emptyenv())
  state$seeds = state$pending = list(); state$anchor = NULL
  state$counts = c(cold = 0L, warm = 0L, fallback = 0L)
  attach_mode_state = function(data, key) {
    cache = data$integration_cache
    cache$mode_state = state; cache$mode_key = key
    for (name in intersect(names(cache), c('selection', 'selection_only', 'with_outcomes', 'outcome')))
      cache[[name]] = attach_mode_state(cache[[name]], paste(key, name, sep = '/'))
    data$integration_cache = cache
    data
  }
  list(prepared = attach_mode_state(prepared, 'joint'), state = state)
}
