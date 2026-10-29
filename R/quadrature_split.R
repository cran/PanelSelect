# Gauss-Legendre with moving cuts for probit and Linear panels with tau fixed at 0.
# The physical (u,v) rule differentiates its finite sum, including moving
# transition cuts, affine nodes and interval lengths. Local index derivatives
# are integrated before expansion through regression designs. No node-by-node
# full parameter tensor, numerical score or integration fallback is used.
#
#   1. shared machinery: panel caches, workspaces, cut partitions and GL axes
#   2. ProbitRE: one probit index (also each factor of the selected panels)
#   3. ProbitRE_ProbitRE: one binary outcome, or the physical (u,v) rule
#   4. ProbitRE_LinearRE: exact outcome marginal, conditioned selection index
#   5. evaluate_split_binary: histories, accuracy gates and diagnostics

# 1. Shared moving-cut machinery ------------------------------------------

split_logsum = function(x) {
  high = max(x)
  if (!is.finite(high)) integration_error('GL quadrature (moving cuts) has no finite mass')
  high + log(sum(exp(x - high)))
}

# Immutable row/design cache; it contains no parameter-dependent quantities.
prepare_split_panels = function(prepared) {
  selected_model = isTRUE(prepared$has_selection)
  active = prepared$active %||% rep(TRUE, prepared$n_rows)
  os = prepared$offset_selection %||% rep(0, prepared$n_rows)
  oy = prepared$offset_outcome %||% rep(0, prepared$n_rows)
  take = function(x, i) {x = x[i, , drop = FALSE]; rownames(x) = NULL; x}
  rows = split(which(active), factor(prepared$group[active], levels = seq_len(prepared$n_panels)))
  panels = lapply(rows, function(i) {
    j = if (selected_model) i[prepared$d[i] == 1] else i
    list(w = if (selected_model) take(prepared$w, i) else NULL,
      x = take(prepared$x, j), selection_sign = if (selected_model) unname(2*prepared$d[i] - 1) else NULL,
      outcome_sign = if (prepared$family == 'probit') unname(2*prepared$y[j] - 1) else NULL,
      outcome = if (prepared$family != 'probit') unname(prepared$y[j]) else NULL,
      log_factorial = if (prepared$family %in% c('poisson', 'pln')) unname(lgamma(prepared$y[j] + 1)) else NULL,
      offset_selection = if (selected_model) unname(os[i]) else NULL,
      offset_outcome = unname(oy[j]))
  })
  if (any(vapply(panels, function(p) any(!is.finite(c(p$w, p$x, p$offset_selection,
      p$offset_outcome, p$outcome, p$log_factorial))) || any(!p$selection_sign %in% c(-1, 1)) ||
      any(!p$outcome_sign %in% c(-1, 1)) || (prepared$family %in% c('poisson', 'pln') &&
        any(p$outcome < 0 | p$outcome != floor(p$outcome))), logical(1))))
    integration_error('GL quadrature (moving cuts) requires finite prepared designs and valid observed outcomes')
  signature = function(x) paste(as.integer(serialize(x, NULL, version = 2)), collapse = ',')
  signatures = vapply(panels, signature, character(1))
  first = which(!duplicated(signatures))
  # Ordered designs remain in these keys: equal indices alone do not imply
  # equal coefficient scores. Never deduplicate a permutation's derivatives.
  keys = lapply(panels[first], function(p) list(
    selection_axis = signature(list(p$w, p$offset_selection)),
    outcome_axis = signature(list(p$x, p$offset_outcome)),
    selection = signature(list(p$w, p$offset_selection, p$selection_sign)),
    outcome = signature(list(p$x, p$offset_outcome, p$outcome_sign, p$outcome, p$log_factorial))))
  fields = c('selection_axis', 'outcome_axis', 'selection', 'outcome')
  ids = lapply(fields, function(name) {
    value = vapply(keys, `[[`, character(1), name)
    as.character(match(value, unique(value)))
  })
  keys = lapply(seq_along(keys), function(k) setNames(lapply(ids, `[[`, k), fields))
  list(panels = panels[first], first = first, map = match(signatures, signatures[first]),
    keys = keys, model_id = prepared$model_id %||%
      resolve_model(prepared$family, prepared$has_selection, prepared$structure)$id,
    family = prepared$family,
    selected_model = selected_model, n_panels = prepared$n_panels)
}

# These small immutable rules depend only on order. Parameter-dependent caches
# below live for ONE evaluation, so a new theta or derivative request cannot
# retrieve stale geometry, products or scores.
split_legendre_rule = local({
  rules = list()
  function(order) {
    key = as.character(order)
    if (is.null(rules[[key]])) {
      rule = statmod::gauss.quad(order, kind = 'legendre')
      if (length(rules) >= 8L) rules <<- tail(rules, 7L)
      rules[[key]] <<- rule
    }
    rules[[key]]
  }
})

split_reuse = function(cache, key, calculate) {
  if (is.null(cache) || is.null(key)) return(calculate())
  if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))
  answer = calculate()
  split_cache_store(cache, key, answer)
  answer
}

# Keep memo storage bounded even when every history/design is different. Large
# individual buffers remain uncached; eviction changes work, never the rule.
split_cache_store = function(cache, key, answer) {
  bytes = as.double(utils::object.size(answer)) + nchar(key, type = 'bytes')
  limit = buffer_bytes
  if (bytes > limit) return(invisible(NULL))
  used = attr(cache, 'bytes') %||% 0
  if (length(cache) >= 64L || used + bytes > limit) {
    rm(list = ls(envir = cache, all.names = TRUE, sorted = FALSE), envir = cache)
    used = 0; attr(cache, 'sizes') = numeric()
  }
  assign(key, answer, envir = cache)
  sizes = attr(cache, 'sizes') %||% numeric(); sizes[key] = bytes
  attr(cache, 'sizes') = sizes
  attr(cache, 'bytes') = used + bytes
  invisible(NULL)
}

split_source_binding = function(prepared) prepared[intersect(c('model_id', 'family', 'has_selection',
  'structure', 'n_rows', 'n_panels', 'group', 'active', 'w', 'x', 'd', 'y',
  'offset_selection', 'offset_outcome'), names(prepared))]

# Call after constructing the evaluator's immutable split_cache. This private
# workspace is captured inside fit_evaluator, not returned in the fitted data.
new_split_workspace = function(prepared) {
  workspace = new.env(parent = emptyenv())
  workspace$source = split_source_binding(prepared)
  workspace$prepared = prepared$split_cache %||% prepare_split_panels(prepared)
  workspace$model_data = prepared$model_data
  workspace
}

split_evaluation_context = function(theta, cache, options, rules, fixed, score, workspace = NULL) {
  create = function() list(partitions = new.env(parent = emptyenv()), axes = new.env(parent = emptyenv()),
    values = new.env(parent = emptyenv()), components = new.env(parent = emptyenv()))
  if (is.null(workspace)) return(create())
  if (!is.environment(workspace)) integration_error('Split evaluator workspace must be private to its evaluator')
  # The prepared binding contains the ordered designs, signs/counts, offsets
  # and model identity. Only this evaluator owns it; there is no global theta
  # memo. Actual rules are included, even when their nominal orders agree.
  settings = list(prepared = cache, options = options, rules = rules,
    fixed = fixed, score = score, parameters = names(theta))
  if (!identical(workspace$settings, settings)) {
    workspace$settings = settings; workspace$context = create()
    workspace$selection = workspace$outcome = workspace$joint = NULL
    workspace$gaussian_key = workspace$gaussian_value = NULL
  }
  selection = theta[c(grep('^selection\\.', names(theta), value = TRUE),
    intersect('sd_re_selection', names(theta)))]
  outcome = theta[c(grep('^outcome\\.', names(theta), value = TRUE),
    intersect(c('sd_re_outcome', 'sd_error'), names(theta)))]
  changed_s = !identical(workspace$selection, selection)
  changed_y = !identical(workspace$outcome, outcome)
  if (changed_s && changed_y) workspace$context = create() else if (
      changed_s || changed_y || !identical(workspace$joint, theta)) {
    pattern = paste0('^(', paste(c(if (changed_s) 'selection', if (changed_y) 'outcome',
      'linear'), collapse = '|'), '):')
    for (memo in workspace$context) {
      remove = grep(pattern, ls(envir = memo, all.names = TRUE, sorted = FALSE), value = TRUE)
      if (length(remove)) {
        rm(list = remove, envir = memo)
        sizes = attr(memo, 'sizes') %||% numeric()
        sizes = sizes[!names(sizes) %in% remove]
        attr(memo, 'sizes') = sizes; attr(memo, 'bytes') = sum(sizes)
      }
    }
  }
  workspace$selection = selection; workspace$outcome = outcome; workspace$joint = theta
  workspace$context
}

# Derivative columns are the local row indices followed by their common SD.
# Exact duplicate cuts are retained as a cluster for the branch-ambiguity
# calculation. Merely nearby cuts are never merged into another finite rule.
split_partition = function(index, scale, bound, score = TRUE) {
  base = sort(unique(c(-bound, seq(-bound, bound, by = 2), bound)))
  n = length(index); cuts = base
  derivatives = if (score) matrix(0, length(base), n + 1L) else NULL
  if (length(index) && scale > 0) {
    levels = c(-4, -2, -1, 0, 1, 2, 4)
    transition = as.vector(outer(-index, levels, '+')/scale)
    keep = which(is.finite(transition) & abs(transition) <= bound)
    cuts = c(cuts, transition[keep])
    if (score) {
      d = matrix(0, length(keep), n + 1L)
      d[cbind(seq_along(keep), (keep - 1L) %% n + 1L)] = -1/scale
      d[, n + 1L] = -transition[keep]/scale
      derivatives = rbind(derivatives, d)
    }
  }
  split_cut_partition(cuts, derivatives, score)
}

# Cut preparation is shared across reference and canonical orders. It preserves
# exact tied cuts and the original near-cut rejection; no interval is merged.
split_cut_partition = function(cuts, derivatives = NULL, score = TRUE) {
  permutation = order(cuts, method = 'radix'); cuts = cuts[permutation]
  if (score) derivatives = derivatives[permutation, , drop = FALSE]
  point = unique(cuts); group = match(cuts, point); width = diff(point)
  if (any(width <= 16*.Machine$double.eps*pmax(1, abs(head(point, -1L)), abs(tail(point, -1L)))))
    integration_error('GL quadrature (moving cuts) has distinct cuts too close for stable interval arithmetic')
  answer = list(points = point, widths = width, clusters = list())
  if (score) {
    clusters = split(seq_along(cuts), factor(group, levels = seq_along(point)))
    representative = t(vapply(clusters, function(j) colMeans(derivatives[j, , drop = FALSE]), numeric(ncol(derivatives))))
    # Domain endpoints never move, including when a transition reaches a bound.
    representative[c(1L, nrow(representative)), ] = 0
    answer$representative = representative
    answer$clusters = lapply(which(lengths(clusters) > 1L), function(k) list(point = k,
      derivatives = derivatives[clusters[[k]], , drop = FALSE], representative = representative[k, ]))
  }
  answer
}

split_axis = function(index, scale, order, bound, score = TRUE, rule = NULL, partition = NULL) {
  if (is.null(rule)) rule = split_legendre_rule(order)
  if (is.null(partition)) partition = split_partition(index, scale, bound, score)
  point = partition$points; width = partition$widths
  segment = rep(seq_along(width), each = order)
  left_fraction = rep((1 - rule$nodes)/2, length(width))
  right_fraction = 1 - left_fraction
  answer = list(nodes = point[segment]*left_fraction + point[segment + 1L]*right_fraction,
    log_weights = log(rep(rule$weights, length(width))) + log(width[segment]/2),
    points = point, widths = width, segment = segment, left_fraction = left_fraction,
    right_fraction = right_fraction, clusters = partition$clusters)
  if (score) {
    representative = partition$representative
    answer$node_derivative = representative[segment, , drop = FALSE]*left_fraction +
      representative[segment + 1L, , drop = FALSE]*right_fraction
    answer$weight_derivative = (representative[segment + 1L, , drop = FALSE] -
      representative[segment, , drop = FALSE])/width[segment]
  }
  answer
}

split_cached_axis = function(index, scale, order, options, score, rule, context = NULL, key = NULL) {
  partition = split_reuse(context$partitions, key,
    function() split_partition(index, scale, options$bound, score))
  split_reuse(context$axes, if (!is.null(key)) paste(key, order, sep = ':'),
    function() split_axis(index, scale, order, options$bound, score, rule, partition))
}

split_expand_local = function(local, design) c(drop(crossprod(design, head(local, -1L))), tail(local, 1L))

# At a tied cut, all possible first-order cut orderings are bounded by the
# residuals of the adjacent finite-rule endpoint derivatives. The boundary
# density here uses the other axis's SAME finite rule. This is a branch bound,
# not a bound on quadrature error or on likelihood score tails.
split_branch_bound = function(axis, mass, flow, boundary_logmass, design, loglik, parameter_map = NULL) {
  answer = rep(0, if (is.null(parameter_map)) ncol(design) + 1L else ncol(parameter_map))
  if (!length(axis$clusters)) return(answer)
  interval_sum = function(x) drop(rowsum(x, axis$segment, reorder = FALSE))
  upper = interval_sum(mass/axis$widths[axis$segment] + flow*axis$right_fraction)
  lower = interval_sum(-mass/axis$widths[axis$segment] + flow*axis$left_fraction)
  for (cluster in axis$clusters) {
    k = cluster$point
    if (is.null(parameter_map)) {
      gradients = cbind(cluster$derivatives[, -ncol(cluster$derivatives), drop = FALSE] %*% design,
        cluster$derivatives[, ncol(cluster$derivatives)])
      center = split_expand_local(cluster$representative, design)
    } else {
      gradients = cluster$derivatives %*% parameter_map
      center = drop(crossprod(parameter_map, cluster$representative))
    }
    radius = apply(abs(sweep(gradients, 2L, center, '-')), 2L, max)
    if (all(radius == 0)) next
    boundary = exp(boundary_logmass(axis$points[k]) - loglik)
    residual = if (k == 1L) abs(lower[1L] + boundary) else if (k == length(axis$points))
      abs(upper[k - 1L] - boundary) else abs(upper[k - 1L] - boundary) + abs(lower[k] + boundary)
    answer = answer + residual*radius
  }
  answer
}

# 2. ProbitRE: one probit index -----------------------------------------

split_probit_values = function(index, sign, scale, nodes, score) {
  z = sweep(outer(nodes, rep(scale, length(index)), '*'), 2L, index, '+')
  z = sweep(z, 2L, sign, '*')
  probability = matrix(pnorm(z, log.p = TRUE), nrow(z), ncol(z))
  result = list(log_product = rowSums(probability))
  if (score) {
    result$partials = sweep(exp(matrix(dnorm(z, log = TRUE), nrow(z), ncol(z)) - probability), 2L, sign, '*')
    result$sum_partials = rowSums(result$partials)
  }
  result
}

# One-dimensional probit panel. score_node is the derivative of log integrand
# with respect to the node, after adding the standard-normal density.
split_probit_panel = function(index, sign, scale, design, order, options, score = TRUE, rule = NULL,
    context = NULL, key = NULL, axis_key = key) {
  component_key = if (!is.null(key)) paste(key, order, sep = ':')
  if (!is.null(context$components) && !is.null(component_key) &&
      exists(component_key, envir = context$components, inherits = FALSE))
    return(get(component_key, envir = context$components, inherits = FALSE))
  if (!length(index)) return(list(loglik = 0, score = rep(0, ncol(design) + 1L),
    ambiguity = rep(0, ncol(design) + 1L), posterior_mean = if (score) 0 else NULL,
    total_nodes = 0L, clusters = 0L, dimension = 0L))
  if (scale == 0) {
    z = sign*index; lp = pnorm(z, log.p = TRUE)
    partials = if (score) sign*exp(dnorm(z, log = TRUE) - lp) else NULL
    return(list(loglik = sum(lp), score = if (score) c(drop(crossprod(design, partials)), 0) else NULL,
      ambiguity = rep(0, ncol(design) + 1L), posterior_mean = if (score) 0 else NULL,
      total_nodes = 0L, clusters = 0L, dimension = 0L))
  }
  axis = split_cached_axis(index, scale, order, options, score, rule, context, axis_key)
  if (length(axis$nodes) > options$max_nodes)
    integration_error(paste0('GL quadrature (moving cuts) node budget exhausted ',
      '(quadrature_options$max_nodes = ', options$max_nodes, ')'))
  values = split_probit_values(index, sign, scale, axis$nodes, score)
  logmass = axis$log_weights + dnorm(axis$nodes, log = TRUE) + values$log_product
  loglik = split_logsum(logmass); gradient = ambiguity = posterior_mean = NULL
  if (score) {
    mass = exp(logmass - loglik)
    flow = mass*(scale*values$sum_partials - axis$nodes)
    local = c(drop(crossprod(mass, values$partials)), sum(mass*axis$nodes*values$sum_partials)) +
      drop(crossprod(axis$node_derivative, flow)) + drop(crossprod(axis$weight_derivative, mass))
    gradient = split_expand_local(local, design)
    boundary = function(z) dnorm(z, log = TRUE) + sum(pnorm(sign*(index + scale*z), log.p = TRUE))
    ambiguity = split_branch_bound(axis, mass, flow, boundary, design, loglik)
    posterior_mean = sum(axis$nodes*mass)
  }
  result = list(loglik = loglik, score = gradient, ambiguity = ambiguity, posterior_mean = posterior_mean,
    total_nodes = length(axis$nodes), clusters = length(axis$clusters), dimension = 1L)
  if (!is.null(context$components) && !is.null(component_key))
    split_cache_store(context$components, component_key, result)
  result
}

# 3. ProbitRE_ProbitRE, tau fixed at 0 ----------------------------------

# For one observed binary outcome and tau explicitly fixed zero, integrate V
# analytically conditional on U. This is a distinct one-dimensional finite rule:
# h=b/d, g=lambda*rho/d, d=sqrt(1+lambda^2*(1-rho^2)). The zero-g integrand
# still has its generally nonzero free-rho derivative through g. Reuse the
# selection cuts and fixed background ONLY: cuts (c-h)/g have a path-dependent
# limit at (h=c,g=0), which would violate the finite-objective contract.
split_single_outcome_panel = function(theta, panel, order, options, score = TRUE,
    clock = function() NULL, rule = NULL, fixed = numeric(), context = NULL) {
  clock()
  if (is.null(panel$w) || nrow(panel$x) != 1L || length(panel$outcome_sign) != 1L ||
      !'cor_error' %in% names(fixed) || fixed[['cor_error']] != 0 || theta[['cor_error']] != 0)
    integration_error('Conditional single-outcome quadrature requires one selected binary outcome and cor_error explicitly fixed at zero')
  selection_names = paste0('selection.', colnames(panel$w))
  outcome_names = paste0('outcome.', colnames(panel$x))
  a = drop(panel$w %*% theta[selection_names]) + panel$offset_selection
  b = drop(panel$x %*% theta[outcome_names]) + panel$offset_outcome
  delta = theta[['sd_re_selection']]; lambda = theta[['sd_re_outcome']]; rho = theta[['cor_re']]
  if (delta < 0 || lambda < 0 || abs(rho) >= 1)
    integration_error('Conditional single-outcome quadrature requires nonnegative scales and interior cor_re')
  for (name in c('sd_re_selection', 'sd_re_outcome')) if (theta[[name]] == 0 && !name %in% names(fixed))
    integration_error('A zero scale must be explicitly fixed for conditional single-outcome quadrature')
  gradient = ambiguity = setNames(rep(0, length(theta)), names(theta))
  finish = function(value, axis = NULL) {
    clock()
    if (!is.finite(value) || (score && any(!is.finite(c(gradient, ambiguity)))))
      integration_error('Conditional single-outcome likelihood or analytical score is not finite')
    list(loglik = value, score = if (score) gradient else NULL,
      ambiguity = if (score) ambiguity else NULL,
      total_nodes = if (is.null(axis)) 0L else length(axis$nodes),
      tied_cut_clusters = if (is.null(axis)) 0L else length(axis$clusters),
      dimension = if (is.null(axis)) 0L else 1L, route = 'conditional_single_outcome',
      tail_relative_bound = if (is.null(axis)) 0 else
        exp(log(2) + pnorm(-options$bound, log.p = TRUE) - value))
  }
  sign_y = panel$outcome_sign
  if (delta == 0) {
    z_s = panel$selection_sign*a; lp_s = pnorm(z_s, log.p = TRUE)
    denominator = sqrt(1 + lambda^2); z_y = sign_y*b/denominator
    lp_y = pnorm(z_y, log.p = TRUE)
    if (score) {
      gradient[selection_names] = drop(crossprod(panel$w,
        panel$selection_sign*exp(dnorm(z_s, log = TRUE) - lp_s)))
      partial = sign_y*exp(dnorm(z_y, log = TRUE) - lp_y)
      gradient[outcome_names] = drop(panel$x)*partial/denominator
      gradient['sd_re_outcome'] = -partial*b*lambda/denominator^3
    }
    return(finish(sum(lp_s) + lp_y))
  }
  denominator = sqrt(1 + lambda^2*(1 - rho^2))
  h = b/denominator; g = lambda*rho/denominator
  if (any(!is.finite(c(denominator, h, g))))
    integration_error('Conditional single-outcome loading is not finite')
  axis = split_cached_axis(a, delta, order, options, score, rule, context,
    if (!is.null(context$keys)) paste0('selection:', context$keys$selection_axis))
  if (length(axis$nodes) > options$max_nodes)
    integration_error(paste0('GL quadrature (moving cuts) node budget exhausted ',
      '(quadrature_options$max_nodes = ', options$max_nodes, ')'))
  s = split_reuse(context$values,
    if (!is.null(context$keys)) paste('selection', context$keys$selection, order, sep = ':'),
    function() split_probit_values(a, panel$selection_sign, delta, axis$nodes, score))
  z_y = sign_y*(h + g*axis$nodes); lp_y = pnorm(z_y, log.p = TRUE)
  logmass = axis$log_weights + dnorm(axis$nodes, log = TRUE) + s$log_product + lp_y
  loglik = split_logsum(logmass)
  if (score) {
    mass = exp(logmass - loglik)
    partial_y = sign_y*exp(dnorm(z_y, log = TRUE) - lp_y)
    flow = mass*(delta*s$sum_partials + g*partial_y - axis$nodes)
    local = c(drop(crossprod(mass, s$partials)), sum(mass*axis$nodes*s$sum_partials),
      sum(mass*partial_y), sum(mass*axis$nodes*partial_y)) +
      c(drop(crossprod(axis$node_derivative, flow)) +
        drop(crossprod(axis$weight_derivative, mass)), 0, 0)
    # Integrate local (a_t,delta,h,g) scores before regression expansion.
    d_denominator = theta*0
    d_denominator['sd_re_outcome'] = lambda*(1 - rho^2)/denominator
    d_denominator['cor_re'] = -lambda^2*rho/denominator
    dh = -h*d_denominator/denominator
    dh[outcome_names] = drop(panel$x)/denominator
    dg = -g*d_denominator/denominator
    dg['sd_re_outcome'] = dg['sd_re_outcome'] + rho/denominator
    dg['cor_re'] = dg['cor_re'] + lambda/denominator
    n = length(a)
    parameter_map = matrix(0, n + 3L, length(theta), dimnames = list(NULL, names(theta)))
    parameter_map[seq_len(n), selection_names] = panel$w
    parameter_map[n + 1L, 'sd_re_selection'] = 1
    parameter_map[n + 2L, ] = dh; parameter_map[n + 3L, ] = dg
    gradient[] = drop(crossprod(parameter_map, local))
    boundary = function(z) dnorm(z, log = TRUE) +
      sum(pnorm(panel$selection_sign*(a + delta*z), log.p = TRUE)) +
      pnorm(sign_y*(h + g*z), log.p = TRUE)
    ambiguity[] = split_branch_bound(axis, mass, flow, boundary, NULL, loglik,
      parameter_map[seq_len(n + 1L), , drop = FALSE])
  }
  finish(loglik, axis)
}

# Selection/outcome factors are calculated once per axis. The generic reduction
# streams finite pairs; regression derivatives are never tensorized.
split_physical_panel = function(theta, panel, order, options, score = TRUE, clock = function() NULL,
    rule = NULL, fixed = numeric(), context = NULL) {
  clock()
  selected_model = !is.null(panel$w)
  if (selected_model && (!'cor_error' %in% names(fixed) || fixed[['cor_error']] != 0 ||
      !'cor_error' %in% names(theta) || theta[['cor_error']] != 0))
    integration_error('Physical GL quadrature (moving cuts) requires cor_error explicitly fixed at zero')
  if (selected_model && isTRUE(options$.single_outcome %||% TRUE) && nrow(panel$x) == 1L &&
      length(panel$outcome_sign) == 1L)
    return(split_single_outcome_panel(theta, panel, order, options, score, clock, rule, fixed, context))
  selection_names = paste0('selection.', colnames(panel$w))
  outcome_names = paste0('outcome.', colnames(panel$x))
  beta = theta[outcome_names]; lambda = theta[['sd_re_outcome']]
  b = drop(panel$x %*% beta) + panel$offset_outcome
  gradient = ambiguity = setNames(rep(0, length(theta)), names(theta))
  finish = function(loglik, nodes, clusters, dimension) {
    if (!is.finite(loglik) || (score && any(!is.finite(c(gradient, ambiguity)))))
      integration_error('GL quadrature (moving cuts) likelihood or analytical score is not finite')
    list(loglik = loglik, score = if (score) gradient else NULL,
      ambiguity = if (score) ambiguity else NULL, total_nodes = nodes,
      tied_cut_clusters = clusters, dimension = dimension, route = 'physical',
      tail_relative_bound = if (dimension) exp(log(2*dimension) +
        pnorm(-options$bound, log.p = TRUE) - loglik) else 0)
  }
  if (!selected_model) {
    result = split_probit_panel(b, panel$outcome_sign, lambda, panel$x, order, options, score, rule,
      context, if (!is.null(context$keys)) paste0('outcome:', context$keys$outcome),
      if (!is.null(context$keys)) paste0('outcome:', context$keys$outcome_axis))
    if (score) {gradient[c(outcome_names, 'sd_re_outcome')] = result$score
      ambiguity[c(outcome_names, 'sd_re_outcome')] = result$ambiguity}
    return(finish(result$loglik, result$total_nodes, result$clusters, result$dimension))
  }
  a = drop(panel$w %*% theta[selection_names]) + panel$offset_selection
  delta = theta[['sd_re_selection']]; rho = theta[['cor_re']]
  # Exact reductions. A free correlation at zero retains its density score.
  if (!length(b) || delta == 0 || lambda == 0 || rho == 0) {
    s = split_probit_panel(a, panel$selection_sign, delta, panel$w, order, options, score, rule,
      context, if (!is.null(context$keys)) paste0('selection:', context$keys$selection),
      if (!is.null(context$keys)) paste0('selection:', context$keys$selection_axis))
    y = split_probit_panel(b, panel$outcome_sign, lambda, panel$x, order, options, score, rule,
      context, if (!is.null(context$keys)) paste0('outcome:', context$keys$outcome),
      if (!is.null(context$keys)) paste0('outcome:', context$keys$outcome_axis))
    if (score) {
      gradient[c(selection_names, 'sd_re_selection')] = s$score
      gradient[c(outcome_names, 'sd_re_outcome')] = y$score
      ambiguity[c(selection_names, 'sd_re_selection')] = s$ambiguity
      ambiguity[c(outcome_names, 'sd_re_outcome')] = y$ambiguity
      if (rho == 0 && delta > 0 && lambda > 0 && length(b)) {
        # This moment is the derivative of the SAME finite physical tensor rule
        # at rho zero. It is generally nonzero even though its value factors.
        gradient['cor_re'] = s$posterior_mean*y$posterior_mean
      }
    }
    return(finish(s$loglik + y$loglik, s$total_nodes + y$total_nodes,
      s$clusters + y$clusters, s$dimension + y$dimension))
  }
  u = split_cached_axis(a, delta, order, options, score, rule, context,
    if (!is.null(context$keys)) paste0('selection:', context$keys$selection_axis))
  v = split_cached_axis(b, lambda, order, options, score, rule, context,
    if (!is.null(context$keys)) paste0('outcome:', context$keys$outcome_axis))
  n = length(u$nodes)*length(v$nodes)
  if (n > options$max_nodes)
    integration_error(paste0('GL quadrature (moving cuts) node budget exhausted ',
      '(quadrature_options$max_nodes = ', options$max_nodes, ')'))
  s = split_reuse(context$values,
    if (!is.null(context$keys)) paste('selection', context$keys$selection, order, sep = ':'),
    function() split_probit_values(a, panel$selection_sign, delta, u$nodes, score))
  y = split_reuse(context$values,
    if (!is.null(context$keys)) paste('outcome', context$keys$outcome, order, sep = ':'),
    function() split_probit_values(b, panel$outcome_sign, lambda, v$nodes, score))
  q = 1 - rho^2
  log_density = function(u, v) -log(2*pi) - log(q)/2 -
    (outer(u^2, rep(1, length(v))) + outer(rep(1, length(u)), v^2) - 2*rho*outer(u, v))/(2*q)
  # The helper sums the SAME finite tensor without materializing its mass.
  # Model density, local scores and branch checks stay in R.
  reduced = bilinear_reduce_cpp(u$log_weights + s$log_product - u$nodes^2/(2*q),
    v$log_weights + y$log_product - v$nodes^2/(2*q), u$nodes, v$nodes, rho/q, score)
  loglik = reduced$value - log(2*pi) - log(q)/2
  if (score) {
    mu = reduced$mu; mv = reduced$mv
    weighted_v = reduced$weighted_v; weighted_u = reduced$weighted_u
    flow_u = delta*s$sum_partials*mu - (u$nodes*mu - rho*weighted_v)/q
    flow_v = lambda*y$sum_partials*mv - (v$nodes*mv - rho*weighted_u)/q
    local_s = c(drop(crossprod(mu, s$partials)), sum(mu*u$nodes*s$sum_partials)) +
      drop(crossprod(u$node_derivative, flow_u)) + drop(crossprod(u$weight_derivative, mu))
    local_y = c(drop(crossprod(mv, y$partials)), sum(mv*v$nodes*y$sum_partials)) +
      drop(crossprod(v$node_derivative, flow_v)) + drop(crossprod(v$weight_derivative, mv))
    gradient[c(selection_names, 'sd_re_selection')] = split_expand_local(local_s, panel$w)
    gradient[c(outcome_names, 'sd_re_outcome')] = split_expand_local(local_y, panel$x)
    gradient['cor_re'] = rho/q + ((1 + rho^2)*sum(u$nodes*weighted_v) -
      rho*(sum(u$nodes^2*mu) + sum(v$nodes^2*mv)))/q^2
    boundary_s = function(z) sum(pnorm(panel$selection_sign*(a + delta*z), log.p = TRUE)) +
      split_logsum(drop(log_density(z, v$nodes)) + v$log_weights + y$log_product)
    boundary_y = function(z) sum(pnorm(panel$outcome_sign*(b + lambda*z), log.p = TRUE)) +
      split_logsum(drop(log_density(u$nodes, z)) + u$log_weights + s$log_product)
    ambiguity[c(selection_names, 'sd_re_selection')] = split_branch_bound(u, mu, flow_u, boundary_s, panel$w, loglik)
    ambiguity[c(outcome_names, 'sd_re_outcome')] = split_branch_bound(v, mv, flow_v, boundary_y, panel$x, loglik)
  }
  clock()
  finish(loglik, n, length(u$clusters) + length(v$clusters), 2L)
}

# 4. ProbitRE_LinearRE, tau fixed at 0 ----------------------------------

split_gaussian_marginal = function(theta, prepared, workspace = NULL) {
  if (is.null(workspace)) return(linear_outcome_marginal(theta, prepared))
  key = list(parameters = theta[c(paste0('outcome.', colnames(prepared$x)),
    prepared$outcome_scale %||% 'sd_re_outcome', 'sd_error')], layout = names(theta))
  if (identical(workspace$gaussian_key, key)) return(workspace$gaussian_value)
  # This evaluator owns its own prepared binding, including changes to y or
  # offsets; do not borrow another integration engine's parameter-only memo.
  prepared$integration_cache = NULL
  value = linear_outcome_marginal(theta, prepared)
  if (as.double(utils::object.size(value)) <= buffer_bytes) {
    workspace$gaussian_key = key; workspace$gaussian_value = value
  } else workspace$gaussian_key = workspace$gaussian_value = NULL
  value
}

# Selected Linear panels, tau fixed zero: integrate the outcome effect exactly.
# U|Y has mean rho*E[V|Y] and variance 1-rho^2+rho^2*Var[V|Y]. The
# remaining one-dimensional split probit rule differentiates this conditioning
# as well as its moving transition cuts; the Gaussian marginal is added once.
split_linear_panel = function(theta, panel, marginal, i, order, options, score, clock, rule, context = NULL) {
  clock()
  alpha_names = paste0('selection.', colnames(panel$w))
  p = length(theta)
  gradient = setNames(marginal$gradient[i, ], names(theta))
  ambiguity = setNames(rep(0, p), names(theta))
  if (is.null(panel$w)) return(list(loglik = marginal$value[i],
    score = if (score) gradient else NULL, ambiguity = if (score) ambiguity else NULL,
    total_nodes = 0L, tied_cut_clusters = 0L, dimension = 0L, tail_relative_bound = 0))
  delta = theta[['sd_re_selection']]; rho = theta[['cor_re']]
  # U given the outcomes: mean rho*mu_v, SD sqrt(1 - rho^2 + rho^2 s_v^2) (ProbitRE_LinearRE.R).
  given = linear_selection_given_outcome(rho, marginal$mean[i], marginal$variance[i],
    marginal$mean_gradient[i, , drop = FALSE], marginal$variance_gradient[i, , drop = FALSE],
    match('cor_re', names(theta)))
  mean_u = given$location; scale_u = given$scale
  a = drop(panel$w %*% theta[alpha_names]) + panel$offset_selection
  design = cbind(panel$w, conditional.location = 1)
  selection = split_probit_panel(a+delta*mean_u, panel$selection_sign,
    delta*scale_u, design, order, options, score, rule, context, paste0('linear:', i))
  if (score) {
    d_mean = given$location_gradient[1L, ]; d_scale = given$scale_gradient[1L, ]
    # The extra design column gives the derivative of a common index shift.
    shift = selection$score[length(alpha_names)+1L]
    scale = tail(selection$score, 1L)
    gradient = gradient + delta*(shift*d_mean + scale*d_scale)
    gradient[alpha_names] = gradient[alpha_names] + selection$score[seq_along(alpha_names)]
    gradient['sd_re_selection'] = gradient['sd_re_selection'] + shift*mean_u + scale*scale_u
    error_shift = selection$ambiguity[length(alpha_names)+1L]
    error_scale = tail(selection$ambiguity, 1L)
    ambiguity = setNames(delta*(error_shift*abs(d_mean) + error_scale*abs(d_scale)), names(theta))
    ambiguity[alpha_names] = ambiguity[alpha_names] + selection$ambiguity[seq_along(alpha_names)]
    ambiguity['sd_re_selection'] = ambiguity['sd_re_selection'] + error_shift*abs(mean_u) + error_scale*scale_u
    names(gradient) = names(ambiguity) = names(theta)
  }
  clock()
  if (!is.finite(marginal$value[i] + selection$loglik) ||
      (score && any(!is.finite(c(gradient, ambiguity)))))
    integration_error('GL quadrature (moving cuts) Linear likelihood or analytical score is not finite')
  list(loglik = marginal$value[i] + selection$loglik,
    score = if (score) gradient else NULL, ambiguity = if (score) ambiguity else NULL,
    total_nodes = selection$total_nodes, tied_cut_clusters = selection$clusters,
    dimension = selection$dimension, route = 'conditional_gaussian', tail_relative_bound = if (selection$dimension)
      exp(log(2) + pnorm(-options$bound, log.p = TRUE) - selection$loglik) else 0)
}

# 5. Driver: histories, accuracy gates and diagnostics ------------------

evaluate_split_binary = function(theta, prepared, control, fixed = numeric(), score = TRUE,
    single_outcome = TRUE) {
  if (!is.logical(score) || length(score) != 1L || is.na(score))
    integration_error('GL quadrature (moving cuts) score must be one logical value')
  options = gl_options(control$nodes, control$quadrature_options %||% list())
  reported_options = options
  # This test switch controls the conditional single-outcome reduction.
  # It is not a public accuracy setting or an automatic numerical fallback.
  if (!is.logical(single_outcome) || length(single_outcome) != 1L || is.na(single_outcome))
    integration_error('Conditional single-outcome dispatch must be one logical value')
  options$.single_outcome = single_outcome
  if (!prepared$family %in% c('probit', 'gaussian') || !identical(prepared$structure, 'panel'))
    integration_error('GL quadrature (moving cuts) currently supports probit and Linear panel models only')
  real_named = function(x, empty = FALSE) is.numeric(x) && !is.complex(x) && all(is.finite(x)) &&
    ((empty && !length(x)) || (!is.null(names(x)) && !anyNA(names(x)) &&
      all(nzchar(names(x))) && !anyDuplicated(names(x))))
  required = parameter_layout(prepared)$name
  if (!real_named(theta) || length(theta) != length(required) || !setequal(names(theta), required) ||
      !real_named(fixed, TRUE) || !all(names(fixed) %in% required) ||
      !all(theta[names(fixed)] == fixed))
    integration_error('GL quadrature (moving cuts) parameters or fixed restrictions do not match the prepared model')
  if (any(theta[intersect(c('sd_re_selection', 'sd_re_outcome'), names(theta))] < 0))
    integration_error('GL quadrature (moving cuts) scales must be nonnegative')
  if (isTRUE(prepared$has_selection)) {
    if (prepared$family == 'probit') split_binary_domain(theta, prepared, fixed, control) else {
      if (!'cor_error' %in% names(fixed) || fixed[['cor_error']] != 0 ||
          theta[['cor_error']] != 0 || abs(theta[['cor_re']]) >= 1)
        integration_error('GL quadrature (moving cuts) for Linear selection requires cor_error explicitly fixed at zero and interior cor_re')
    }
  }
  if ('sd_error' %in% names(theta) && theta[['sd_error']] <= 0)
    integration_error('GL quadrature (moving cuts) for Linear requires a positive residual scale')
  for (name in intersect(c('sd_re_selection', 'sd_re_outcome'), names(theta)))
    if (theta[[name]] == 0 && !name %in% names(fixed))
      integration_error('A zero scale must be explicitly fixed for GL quadrature (moving cuts)')
  workspace = prepared$split_workspace
  if (!is.null(workspace)) {
    if (!is.environment(workspace)) integration_error('Split evaluator workspace must be private to its evaluator')
    source = split_source_binding(prepared)
    if (!identical(workspace$source, source)) {
      prepared$model_data = prepared$integration_cache = NULL
      prepared = prepare_model_data(prepared)
      workspace$source = source; workspace$prepared = prepare_split_panels(prepared)
      workspace$model_data = prepared$model_data; workspace$score_floors = NULL
    }
    if (is.null(workspace$model_data)) workspace$model_data = prepare_model_data(prepared)$model_data
    prepared$model_data = workspace$model_data
    cache = workspace$prepared
  } else cache = prepared$split_cache %||% prepare_split_panels(prepared)
  orders = c(options$reference_nodes, as.integer(control$nodes))
  rules = lapply(orders, split_legendre_rule)
  floors = score_floors(prepared, theta, workspace)
  context = split_evaluation_context(theta, cache, options, rules, fixed, score, workspace)
  marginal = if (prepared$family == 'gaussian') split_gaussian_marginal(theta, prepared, workspace) else NULL
  if (!is.null(marginal) && (any(!is.finite(marginal$value)) ||
      (score && any(!is.finite(marginal$gradient)))))
    integration_error('GL quadrature (moving cuts) Linear Gaussian marginal or analytical score is not finite')
  values = numeric(length(cache$panels)); changes = values
  gradients = matrix(0, length(cache$panels), length(theta), dimnames = list(NULL, names(theta)))
  reference = error = gradients; diagnostics = vector('list', length(cache$panels))
  started = proc.time()[['elapsed']]
  fail = function(message, details = NULL) stop(structure(list(message = paste('GL quadrature (moving cuts):', message),
    call = NULL, diagnostics = details), class = c('panelselect_split_accuracy_error',
      'panelselect_integration_error', 'error', 'condition')))
  for (k in seq_along(cache$panels)) {
    context$keys = cache$keys[[k]]
    tick = proc.time()[['elapsed']]
    clock = function() if (proc.time()[['elapsed']] - tick > options$seconds_per_panel ||
      proc.time()[['elapsed']] - started > options$seconds_total) fail('time budget exhausted (quadrature_options$seconds_per_panel, $seconds_total)', list(history = k))
    evaluate = function(j) if (is.null(marginal))
      split_physical_panel(theta, cache$panels[[k]], orders[j], options, score, clock, rules[[j]], fixed, context) else
      split_linear_panel(theta, cache$panels[[k]], marginal, cache$first[k], orders[j], options, score, clock, rules[[j]], context)
    predecessor = evaluate(1L); canonical = evaluate(2L)
    changes[k] = abs(canonical$loglik - predecessor$loglik)
    if (any(!is.finite(c(canonical$tail_relative_bound, predecessor$tail_relative_bound))) ||
        max(canonical$tail_relative_bound, predecessor$tail_relative_bound) > options$tail_tol)
      fail('probability tail bound unresolved (quadrature_options$tail_tol, $bound)', list(history = k, canonical = canonical, predecessor = predecessor))
    if (changes[k] > options$panel_loglik_tol)
      fail('per-panel likelihood refinement unresolved (quadrature_options$panel_loglik_tol)', list(history = k, change = changes[k]))
    if (score) {
      canonical$score[names(fixed)] = predecessor$score[names(fixed)] = 0
      canonical$ambiguity[names(fixed)] = predecessor$ambiguity[names(fixed)] = 0
      error[k, ] = abs(canonical$score - predecessor$score) + canonical$ambiguity + predecessor$ambiguity
      if (any(error[k, ] > options$panel_score_tol*(floors + abs(canonical$score))))
        fail('per-panel analytical score refinement or cut-branch ambiguity unresolved (quadrature_options$panel_score_tol)',
          list(history = k, error = error[k, ], canonical = canonical, predecessor = predecessor))
      gradients[k, ] = canonical$score; reference[k, ] = predecessor$score
    }
    values[k] = canonical$loglik
    diagnostics[[k]] = list(canonical_order = orders[2], reference_order = orders[1],
      dimension = canonical$dimension, nodes = canonical$total_nodes, reference_nodes = predecessor$total_nodes,
      route = canonical$route %||% 'physical',
      loglik_change = changes[k], maximum_score_error = if (score) max(c(0, error[k, ])) else NA_real_,
      maximum_branch_ambiguity = if (score) max(c(0, canonical$ambiguity, predecessor$ambiguity)) else NA_real_,
      tied_cut_clusters = canonical$tied_cut_clusters, tail_relative_bound = canonical$tail_relative_bound,
      elapsed_seconds = proc.time()[['elapsed']] - tick)
  }
  absolute_change = sum(changes[cache$map])
  natural_error = colSums(error[cache$map, , drop = FALSE])
  jacobian = setNames(rep(1, length(theta)), names(theta))
  sd = intersect(c('sd_re_selection', 'sd_re_outcome', 'sd_error'), names(theta)); jacobian[sd] = theta[sd]
  correlations = intersect(c('cor_re', 'cor_error'), names(theta)); jacobian[correlations] = 1 - theta[correlations]^2
  jacobian[names(fixed)] = 0
  # Optimizer-coordinate errors; each budget scales with its score floor.
  aggregate_error = natural_error*abs(jacobian)
  score_error = max(c(0, aggregate_error))
  if (absolute_change > options$total_loglik_tol ||
      (score && any(aggregate_error > options$total_score_tol*floors)))
    fail('aggregate likelihood or score refinement budget unresolved (quadrature_options$total_loglik_tol, $total_score_tol)', list(loglik = absolute_change, score = score_error))
  routes = vapply(diagnostics, function(x) x$route, character(1))
  list(loglik = values[cache$map], score = if (score) gradients[cache$map, , drop = FALSE] else NULL,
    reference_score = if (score) reference[cache$map, , drop = FALSE] else NULL,
    diagnostics = list(method = 'split_legendre', status = 'PASS', options = reported_options,
      coordinates = if (!is.null(marginal)) 'exact Gaussian outcome; standardized conditioned selection effect' else
        if (all(routes == 'conditional_single_outcome')) 'standardized selection effect; single outcome integrated exactly' else
          if (any(routes == 'conditional_single_outcome')) 'physical and conditional single-outcome paths' else 'physical',
      canonical_order = orders[2], reference_order = orders[1],
      nodes = vapply(diagnostics, function(x) x$nodes, numeric(1)),
      dimension = vapply(diagnostics, function(x) x$dimension, integer(1)),
      routes = routes,
      unique_histories = length(cache$panels), panel_to_history = cache$map, histories = diagnostics,
      score_checked = isTRUE(score), check_scope = if (score)
        'likelihood/refinement/tail/analytic score' else 'likelihood/refinement/tail only',
      aggregate_absolute_loglik_change = absolute_change, aggregate_absolute_score_change = score_error,
      aggregate_score_units = 'optimizer coordinates', gradient_method = 'analytic',
      score_basis = 'Derivative of the finite split rule including cuts, nodes and weights; at exact crossings a symmetric analytical branch score with an explicit ambiguity budget.',
      fallback = 'none', probability_floor = 'none',
      accuracy_limitation = 'Order refinement is not a rigorous remainder certificate. Probability tails do not bound score tails. Physical Gaussian ridges near absolute rho one may fail refinement.'))
}

# Selected binary panels with cor_error explicitly fixed at zero.
split_binary_domain = function(theta, prepared, fixed, control) {
  fail = function(message) integration_error(paste('Split GL binary integration:', message))
  if (!identical(prepared$family, 'probit') || !isTRUE(prepared$has_selection) ||
      !identical(prepared$structure, 'panel'))
    fail('only selected binary panel models are supported')
  if (!'cor_error' %in% names(fixed) || fixed[['cor_error']] != 0)
    fail('cor_error must be explicitly fixed at zero; a free value currently at zero is unsupported')
  real_named = function(x, empty = FALSE) is.numeric(x) && !is.complex(x) && all(is.finite(x)) &&
    ((empty && !length(x)) || (!is.null(names(x)) && !anyNA(names(x)) &&
       all(nzchar(names(x))) && !anyDuplicated(names(x))))
  required = c(paste0('selection.', colnames(prepared$w)),
    paste0('outcome.', colnames(prepared$x)), 'sd_re_selection', 'sd_re_outcome', 'cor_re', 'cor_error')
  if (!real_named(theta) || length(theta) != length(required) || !setequal(names(theta), required) ||
      !real_named(fixed, TRUE) || !all(names(fixed) %in% required) ||
      !all(theta[names(fixed)] == fixed)) fail('natural parameters or fixed restrictions do not match the design')
  if (theta[['cor_error']] != 0 || any(theta[c('sd_re_selection', 'sd_re_outcome')] < 0) ||
      abs(theta[['cor_re']]) >= 1) fail('parameters are outside the supported domain')
  invisible(TRUE)
}
