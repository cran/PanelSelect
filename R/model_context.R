# Inputs to the numerical kernels.
#
# A prepared data set owns observations, designs, panel membership and fixed
# model identity. prepare_model_data caches its immutable row blocks. A model
# card supplies the functions for that identity; make_model_context binds them
# to one natural parameter vector theta, producing the point m used in model
# files. It contains indices a, b, score columns and affine person-effect maps.
#
# Read model_dimension before the maps: restrictions determine the integration
# dimension, and the Linear point hook can condition the maps on observed
# outcomes. subset_panels and selection_component construct derived data sets
# while preserving parameter identities and each panel's row order.

# A model card binds metadata to model-specific or shared functions. likelihood
# defaults to integrate_person_effects; standalone cards use integrate_standalone,
# and selection models can supply entries with additional reductions. panel is
# the explicit kernel at person-effect nodes; outcome_rows (standalone models
# and the nested fixed rule) may reuse another model's row function and gives the
# outcome rows' log likelihood and partials at given indices B; point may
# replace the standard person-effect maps (Linear conditions them on the
# observed outcomes); exact_outcome gives an exact outcome factor for the
# independent factorization (Linear); predictions gives the population, joint
# and selected means.
# These runtime lists are private; fitted objects store serializable model IDs.
model_definition = function(id, file, structure, family, has_selection, parameters,
                            panel, predictions, likelihood = integrate_person_effects,
                            outcome_rows = NULL, point = NULL, exact_outcome = NULL,
                            positive_scales = character(), zero_correlations = list(),
                            validate_design = NULL) {
  list(id = id, file = file, structure = structure, family = family,
       has_selection = has_selection, parameters = parameters,
       likelihood = likelihood, panel = panel, outcome_rows = outcome_rows, point = point,
       exact_outcome = exact_outcome, predictions = predictions,
       positive_scales = positive_scales, zero_correlations = zero_correlations,
       validate_design = validate_design,
       capabilities = list(history = structure == 'panel',
         random_effects = structure == 'panel', inner = family == 'pln'))
}

resolve_model = function(family, has_selection, structure = 'panel', model_id = NULL) {
  if (identical(model_id, 'probit_poisson') ||
      (is.null(model_id) && identical(structure, 'cross_sectional') &&
       isTRUE(has_selection) && identical(family, 'poisson')))
    stop('family = "poisson" requires a panel ID; for cross-sectional count selection, use family = "count" or "pln"',
      call. = FALSE)
  if (is.null(model_id)) {
    # Structure is explicit metadata; its absence never selects a panel model.
    if (is.null(structure))
      stop('Prepared data have missing metadata: neither a model identity nor a structure', call. = FALSE)
    if (structure == 'cross_sectional') {
      if (!has_selection) stop('Cross-sectional models require selection', call. = FALSE)
      model_id = switch(family, probit = 'probit_probit', gaussian = 'probit_linear',
                        pln = 'probit_pln')
    } else if (has_selection) {
      model_id = switch(family, probit = 'probit_re_probit_re', gaussian = 'probit_re_linear_re',
                        poisson = 'probit_re_poisson_re', pln = 'probit_re_pln_re')
    } else model_id = switch(family, probit = 'probit_re', poisson = 'poisson_re', pln = 'pln_re')
  }
  if (is.null(model_id)) stop('Unsupported model identity', call. = FALSE)
  model = switch(model_id,
    probit_re = probit_re_model(),
    poisson_re = poisson_re_model(),
    pln_re = pln_re_model(),
    probit_re_probit_re = probit_re_probit_re_model(),
    probit_re_linear_re = probit_re_linear_re_model(),
    probit_re_poisson_re = probit_re_poisson_re_model(),
    probit_re_pln_re = probit_re_pln_re_model(),
    probit_probit = probit_probit_model(),
    probit_linear = probit_linear_model(),
    probit_pln = probit_pln_model())
  if (is.null(model)) stop('Unsupported model identity', call. = FALSE)
  model
}

# Restrictions, dimensions and the model point ---------------------------------

# An absent parameter is structurally zero. A present parameter counts as
# fixed zero only through fixed, not merely because theta currently equals 0.
is_fixed_zero = function(name, theta, fixed) {
  !name %in% names(theta) || (name %in% names(fixed) && fixed[[name]] == 0)
}

`%||%` = function(x, y) if (is.null(x)) y else x

# q counts coordinates used by the likelihood, not always both model effects.
# Positive person-effect scales activate their coordinates; exact Linear
# conditioning with tau fixed at 0 integrates only the selection coordinate.
model_dimension = function(theta, prepared, control, fixed) {
  if (identical(prepared$structure, 'cross_sectional')) return(0L)
  outcome = theta[[prepared$outcome_scale %||% 'sd_re_outcome']] > 0
  if (!prepared$has_selection) return(as.integer(outcome))
  selection = theta[['sd_re_selection']] > 0
  if (prepared$family == 'gaussian' && isTRUE(control$exact_reductions) &&
      is_fixed_zero('cor_error', theta, fixed)) return(as.integer(selection))
  as.integer(selection) + as.integer(outcome)
}

# The one remaining latent effect when model_dimension() is 1. Length-two
# nodes give the selection effect's order first and the outcome effect's second.
single_effect = function(theta, prepared) {
  if (!prepared$has_selection)
    return(if (identical(prepared$outcome_scale, 'sd_re_selection')) 'selection' else 'outcome')
  if (theta[['sd_re_selection']] > 0) 'selection' else 'outcome'
}

# Cache row blocks once per prepared data set. ptr contains zero-based group
# boundaries for reductions; rows retain their original within-panel order.
# The selected block alone supplies outcomes to selection-model row kernels.
prepare_model_data = function(prepared) {
  if (!is.null(prepared$model_data)) return(prepared)
  if (is.null(prepared$model_id)) prepared$model_id =
    resolve_model(prepared$family, prepared$has_selection, prepared$structure)$id
  block = function(rows) {
    count = tabulate(prepared$group[rows], prepared$n_panels)
    list(rows = rows, group = prepared$group[rows], ptr = c(0L, cumsum(count)),
         x = prepared$x[rows, , drop = FALSE], w = prepared$w[rows, , drop = FALSE],
         y = prepared$y[rows], count = count)
  }
  active = prepared$active %||% rep(TRUE, prepared$n_rows)
  selected = block(which(prepared$d == 1 & active))
  selected$sum_x = group_sum(selected$x, selected$ptr)
  prepared$model_data = list(selected = selected,
    unselected = block(which(prepared$d == 0 & active)),
    active = block(which(active)),
    outcome_names = if (ncol(prepared$x)) paste0(prepared$outcome_prefix %||% 'outcome.', colnames(prepared$x)) else character(),
    selection_names = if (ncol(prepared$w)) paste0('selection.', colnames(prepared$w)) else character(),
    log_factorial = if (prepared$family %in% c('poisson', 'pln')) lgamma(prepared$y+1) else NULL)
  prepared
}

# Bind one natural parameter point. Score columns here index theta; the kernel
# binder adds a latent-column offset when it requests mixed derivatives.
# indices_only is the smaller context for collapsed standalone fixed rules.
# The optional model point hook supplies conditioned maps and an exact marginal;
# otherwise standard_person_effects supplies the prior maps.
make_model_context = function(theta, prepared, control, fixed, indices_only = FALSE) {
  prepared = prepare_model_data(prepared)
  definition = resolve_model(prepared$family, prepared$has_selection,
    prepared$structure, prepared$model_id)
  # Natural columns of the local structural coordinates: error is sd_error and
  # tau is cor_error in every model.
  structural_columns = c(error = match('sd_error', names(theta)),
    tau = match('cor_error', names(theta)))
  data = prepared$model_data; n = prepared$n_panels; p = length(theta)
  q = model_dimension(theta, prepared, control, fixed)
  parameter = function(name) if (name %in% names(theta)) unname(theta[[name]]) else 0
  ix = match(data$outcome_names, names(theta)); iw = match(data$selection_names, names(theta))
  predictor = function(x, index, offset) if (length(index)) drop(x %*% theta[index])+offset else offset
  b = predictor(prepared$x, ix, prepared$offset_outcome)
  a = if (prepared$has_selection) predictor(prepared$w, iw, prepared$offset_selection) else NULL
  delta = parameter('sd_re_selection'); lambda_name = prepared$outcome_scale %||% 'sd_re_outcome'
  lambda = parameter(lambda_name); error = parameter('sd_error'); rho = parameter('cor_re')
  tau = parameter('cor_error'); tau_zero = is_fixed_zero('cor_error', theta, fixed)
  gaussian = prepared$family == 'gaussian' && definition$structure == 'panel' &&
    isTRUE(control$exact_reductions)
  il = match(lambda_name, names(theta)); ir = match('cor_re', names(theta))
  # Natural parameter columns are resolved once, including selection-only aliases.
  column = function(name) which(names(theta) == name)
  score_columns = list(alpha = iw, beta = ix, delta = column('sd_re_selection'),
    lambda = column(lambda_name), rho = column('cor_re'), tau = column('cor_error'),
    sigma = if (prepared$family == 'gaussian') column('sd_error') else integer(),
    gamma = if (prepared$family == 'pln') column('sd_error') else integer())
  # A PLN observation effect with a positive scale needs its observation rule
  # (quadrature_inner.R); a fixed gamma = 0 removes it.
  inner = if (definition$capabilities$inner && error > 0) observation_rule(theta, prepared,
    control, list(a = a, b = b, delta = delta, lambda = lambda, lambda_name = lambda_name,
      rho = rho, alpha = iw, beta = ix))
  measure_recorder = control$.staged_recorder
  measure_component = if (!is.null(measure_recorder) &&
      isTRUE(measure_recorder$has_selection) && !prepared$has_selection) {
    if (identical(prepared$outcome_scale, 'sd_re_selection')) 'selection' else 'outcome'
  } else 'joint'
  # Fixed standalone integration needs indices, not latent effect maps.
  if (indices_only) return(list(data = data, definition = definition,
    structural_columns = structural_columns, score_columns = score_columns, q = q, b = b, ix = ix,
    lambda = lambda, lambda_index = il, error = error, inner = inner, gaussian = FALSE,
    staged_recorder = measure_recorder, staged_reference = isTRUE(control$.staged_reference),
    staged_component = measure_component))
  # Person effects are affine in z: location + loading*z per unit, with their
  # parameter derivatives (q columns of z, p natural parameters).
  effect = function() list(location = numeric(n), loading = matrix(0, n, q),
    location_gradient = matrix(0, n, p), loading_gradient = array(0, c(n, q, p)))
  m = list(data = data, definition = definition,
    structural_columns = structural_columns, score_columns = score_columns, q = q, b = b, a = a,
       d = prepared$d, u = effect(), v = effect(), marginal = NULL,
       ix = ix, iw = iw, delta = delta, lambda = lambda, lambda_index = il,
       error = error, tau = tau, tau_zero = tau_zero, gaussian = gaussian, inner = inner,
       staged_recorder = measure_recorder, staged_reference = isTRUE(control$.staged_reference),
       staged_component = measure_component, rho = rho, rho_index = ir,
       has_u = prepared$has_selection && delta > 0, has_v = lambda > 0)
  # The standard maps, or the model's own (Linear conditions on the outcomes).
  m = (definition$point %||% standard_person_effects)(m, theta, prepared)
  # Most prior maps depend on no parameter, or only the person correlation.
  # Retain exact nonzero derivative columns when a first-order kernel expands z.
  active_columns = function(map) {
    which(colSums(map$location_gradient != 0 | is.na(map$location_gradient)) +
      colSums(matrix(map$loading_gradient != 0 | is.na(map$loading_gradient), n*q, p)) > 0)
  }
  m$u$columns = active_columns(m$u); m$v$columns = active_columns(m$v)
  m
}

# Independent standard-normal z: u = z1 and v = rho*z1 + sqrt(1-rho^2)*z2 when
# both effects are active (a single active effect is z1), with dv/drho.
standard_person_effects = function(m, theta, prepared) {
  if (m$has_u) m$u$loading[, 1L] = 1
  if (m$has_v) {
    if (m$has_u) {
      m$v$loading[, 1L] = m$rho; m$v$loading[, 2L] = sqrt(1-m$rho^2)
      if (!is.na(m$rho_index)) {
        m$v$loading_gradient[, 1L, m$rho_index] = 1
        m$v$loading_gradient[, 2L, m$rho_index] = -m$rho/sqrt(1-m$rho^2)
      }
    } else m$v$loading[, 1L] = 1
  }
  m
}

subset_panels = function(prepared, panels) {
  if (is.null(panels)) panels = integer()
  if (!is.numeric(panels) || is.complex(panels) || !is.null(dim(panels)) ||
      anyNA(panels) || any(!is.finite(panels)) || any(panels != trunc(panels)) ||
      any(panels < 1 | panels > prepared$n_panels) || anyDuplicated(panels))
    stop('panels must be distinct finite integer indices within the prepared panels',
         call. = FALSE)
  panels = as.integer(panels)
  rows = which(prepared$group %in% panels)
  # Group pointers describe contiguous blocks in the requested panel order.
  # Preserve each panel's row order, including unequal and empty panel lengths.
  if (is.unsorted(panels))
    rows = rows[order(match(prepared$group[rows], panels), method = 'radix')]
  ans = prepared
  for (name in c('x', 'w')) ans[[name]] = prepared[[name]][rows, , drop = FALSE]
  for (name in c('y', 'd', 'offset_outcome', 'offset_selection', 'selected', 'active',
                 '.staged_rows'))
    if (!is.null(prepared[[name]])) ans[[name]] = prepared[[name]][rows]
  ans$group = match(prepared$group[rows], panels)
  ans$n_rows = length(rows); ans$n_panels = length(panels)
  ans$n_units = ans$n_panels
  ans$group_ptr = c(0L, cumsum(tabulate(ans$group, nbins = ans$n_panels)))
  ans$id_values = prepared$id_values[panels]
  if (!is.null(prepared$.staged_panel_ids)) ans$.staged_panel_ids = prepared$.staged_panel_ids[panels]
  if (!is.null(prepared$original_index)) {
    ans$original_index = prepared$original_index[rows]
    ans$keep = sort(ans$original_index)
    ans$inverse_index = order(ans$original_index)
  }
  ans$model_data = ans$integration_cache = NULL
  ans
}

selection_component = function(prepared) {
  ans = prepared
  ans$model_id = 'probit_re'
  ans$family = 'probit'; ans$has_selection = FALSE
  ans$x = prepared$w; ans$y = prepared$d; ans$d = rep(1, prepared$n_rows)
  ans$offset_outcome = prepared$offset_selection
  ans$outcome_prefix = 'selection.'; ans$outcome_scale = 'sd_re_selection'
  ans$active = rep(TRUE, prepared$n_rows)
  ans$model_data = ans$integration_cache = NULL
  ans
}
