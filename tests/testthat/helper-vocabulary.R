# Historical vocabulary constructed from current fits for tests of the
# normalizers (test-saved-vocabulary.R, test-normalizers.R). Actual historical
# objects are in fixtures/saved-vocabulary.R (test-saved-objects.R).

# Rewrite a current fit in the historical schema: rule values adaptive, fixed and split;
# the inner_* control fields; one option list per rule with the earlier option
# names; reduction and stationarity_tol; the fit settings among the control
# fields; sd_selection, sd_outcome and cor_id; loglik_by_panel; meshes labelled split
# with the earlier binding; and a call record in the earlier argument order.
old_vocabulary = function(fit, formula_record = FALSE) {
  rule = c(agh = 'adaptive', gh = 'fixed', gl = 'split', gk = 'gk')
  option_names = c(panel_loglik_tol = 'loglik_tol', panel_score_tol = 'score_tol',
    total_loglik_tol = 'aggregate_loglik_tol', total_score_tol = 'aggregate_score_tol',
    max_cells = 'max_intervals')
  old_options = function(x) {
    hit = names(x) %in% names(option_names)
    names(x)[hit] = option_names[names(x)[hit]]
    x
  }
  per_rule = function(options, current, wanted) if (identical(current, wanted)) old_options(options) else list()
  parameters = c(sd_re_selection = 'sd_selection', sd_re_outcome = 'sd_outcome', cor_re = 'cor_id')
  old_names = function(x) {
    hit = x %in% names(parameters)
    x[hit] = parameters[x[hit]]
    x
  }
  rename = function(x) {
    if (is.data.frame(x)) {
      for (j in seq_along(x)) if (is.character(x[[j]])) x[[j]] = old_names(x[[j]])
      return(x)
    }
    if (is.list(x)) {
      for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] = rename(x[[i]])
    } else if (is.character(x)) x = old_names(x)
    if (!is.null(names(x))) names(x) = old_names(names(x))
    if (!is.null(dim(x)) && !is.null(dimnames(x)))
      dimnames(x) = lapply(dimnames(x), function(d) if (is.null(d)) d else old_names(d))
    x
  }
  # Meshes labelled their GL rule 'split'.
  relabel = function(x) {
    if (!is.list(x)) return(x)
    if (identical(x$method, 'gl')) x$method = 'split'
    for (i in seq_along(x)) if (is.list(x[[i]])) x[[i]] = relabel(x[[i]])
    x
  }
  control = fit$control
  outer = control$quadrature; inner = control$observation_quadrature
  earlier = list(quadrature = rule[[outer]], nodes = control$nodes,
    inner_nodes = control$observation_nodes, optimizer = control$optimizer,
    method = control$method, maxit = control$maxit, reltol = control$reltol,
    factr = control$factr, mode_maxit = control$mode_maxit, mode_tol = control$mode_tol,
    vcov = control$vcov, trace = control$trace, na_action = control$na_action,
    init = control$init, fixed_geometry = NULL,
    reduction = control$exact_reductions, inner_quadrature = rule[[inner]],
    binary_integration = 'gh', checked = list(),
    split = per_rule(control$quadrature_options, outer, 'gl'),
    gk = per_rule(control$quadrature_options, outer, 'gk'),
    inner_split = per_rule(control$observation_options, inner, 'gl'),
    inner_gk = per_rule(control$observation_options, inner, 'gk'),
    parscale = control$parscale, stationarity_tol = control$stationarity_warning)
  mesh = control$integration_mesh
  if (!is.null(mesh)) {
    b = mesh$binding
    head = b[seq_len(match('fixed', names(b)) - 1L)]
    mesh$binding = c(head, list(fixed = b$fixed, reduction = b$exact_reductions,
      outer = rule[[b$outer]], nodes = b$nodes, inner = rule[[b$inner]],
      physical_pln = b$physical_pln, inner_nodes = b$observation_nodes,
      split = per_rule(b$quadrature_options, b$outer, 'gl'),
      gk = per_rule(b$quadrature_options, b$outer, 'gk'),
      inner_split = per_rule(b$observation_options, b$inner, 'gl'),
      inner_gk = per_rule(b$observation_options, b$inner, 'gk')))
    earlier$integration_mesh = relabel(mesh)
  }
  attributes(earlier) = c(list(names = names(earlier)),
    attributes(control)[setdiff(names(attributes(control)), 'names')])
  old = fit
  old$control = earlier
  names(old)[names(old) == 'loglik_by_unit'] = 'loglik_by_panel'
  for (name in setdiff(names(old), c('data', 'prepared', 'call')))
    if (!is.null(old[[name]])) old[[name]] = rename(old[[name]])
  # The historical call record: arguments in its formal order, settings in the control.
  call = as.list(fit$call)
  keep = intersect(c('outcome', 'data', 'id', 'family', 'selection', 'start', 'fixed', 'control',
    'subset'), names(call))
  call = c(call[1L], call[keep])
  if (formula_record) names(call)[names(call) == 'outcome'] = 'formula'
  old$call = as.call(call)
  old
}
