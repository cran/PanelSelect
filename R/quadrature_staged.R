# Staged Gauss-Legendre integration with positive GL meshes.
# Partitions are prepared between optimization stages; every value, score and
# information stencil within a stage uses the same mesh. Model files supply
# observation likelihoods and natural first scores; shared R primitives and
# assembly supply their higher mixed derivatives. Outer mesh evaluation uses
# make_panel_kernel, including its standalone first-score fast path. Eligible
# selected PLN models can use the separate physical (u,v) evaluator below.
#   quadrature_staged.R          routing, bindings, cells, tails, outer evaluation
#   quadrature_staged_inner.R    PLN observation meshes, records and probes
#   quadrature_staged_pln.R      physical (u,v) rule for selected PLN, tau fixed 0
#   quadrature_staged_prepare.R  seeding and refining the next stage's mesh
# Lifecycle: prepare_staged_integration seeds/refines a mesh between stages;
# evaluate_staged_integration reads it at one parameter point. A fresh recorder
# accumulates diagnostics for that evaluation; a bounded evaluator workspace
# reuses only parameter-independent cell rules. Saved meshes bind to data,
# model, controls and fixed restrictions, not to a live recorder or workspace.
# Canonical/reference differences and tail bounds diagnose this finite mesh;
# they do not turn its score into the exact integral's score.

# 1. Routing, controls and failure conditions ----------------------------

uses_staged_integration = function(prepared, control, fixed = numeric()) {
  if (isTRUE(control$.staged_evaluation)) return(FALSE)
  if (identical(prepared$family, 'pln') &&
      (control$observation_quadrature %||% 'gh') == 'gl' &&
      !('sd_error' %in% names(fixed) && fixed[['sd_error']] == 0)) return(TRUE)
  if (!uses_gl_quadrature(control)) return(FALSE)
  # Preserve the existing physical moving-cut panel routes by restrictions,
  # including free rho evaluated at zero. A free tau at zero stays coupled.
  # 'old' names the supported moving-cut route, not an obsolete implementation.
  old = identical(prepared$structure, 'panel') &&
    prepared$family %in% c('probit', 'gaussian') &&
    (!prepared$has_selection || ('cor_error' %in% names(fixed) && fixed[['cor_error']] == 0))
  !old
}

staged_outer_method = function(control)
  if (uses_gl_quadrature(control)) 'gl' else control$quadrature

# PLN observation meshes are staged only when the observation error is active.
# 'new_inner' denotes this rule choice, not whether a mesh has been initialized.
staged_new_inner = function(prepared, control, theta)
  prepared$family == 'pln' && theta[['sd_error']] > 0 &&
    (control$observation_quadrature %||% 'gh') == 'gl'

# Native GH has no public outer time controls. Preserve the finite overall
# default of the GL rule independently of the user's active observation-work
# limits.
staged_total_seconds = function(control, options)
  if (control$quadrature %in% c('gh', 'agh') && !uses_gl_quadrature(control))
    max(gl_options(control$observation_nodes)$seconds_total,
      staged_options(control, TRUE)$seconds_total) else options$seconds_total

# A time guard that first checks its enclosing guard.
staged_clock = function(seconds, message, parent = NULL) {
  started = proc.time()[['elapsed']]
  function() {
    if (!is.null(parent)) parent()
    if (proc.time()[['elapsed']]-started > seconds) integration_error(message)
    invisible(NULL)
  }
}

# GL keeps its bound: the domain is never widened. An omitted tail that the
# fixed domain leaves unresolved stops the preparation with the option that
# sets the bound, `name` (staged_options_name(bound = TRUE)).
staged_bound_failure = function(message, bound, name, theta, diagnostics)
  staged_mesh_required(paste0(message, ' (maximum bound ', format(bound), '; raise ', name, ')'),
    theta, diagnostics)

# The control option list that sets a staged limit: quadrature_options for an
# outer GL rule, otherwise observation_options (an outer GH rule with a staged
# observation mesh uses the observation options, as staged_options() does).
# With `bound = TRUE`, the option that sets the domain bound.
staged_options_name = function(control, inner = FALSE, bound = FALSE) {
  if (!inner && staged_outer_method(control) != 'gl') return(staged_options_name(control, TRUE, bound))
  paste0(if (inner) 'observation_options' else 'quadrature_options', if (bound) '$bound')
}

staged_options = function(control, inner = FALSE) {
  method = if (inner) control$observation_quadrature else staged_outer_method(control)
  # Rule names are checked here, so no name falls through to another rule.
  if (!isTRUE(method %in% c('agh', 'gh', 'gl')) || (inner && method != 'gl'))
    stop('Internal error: no staged ', if (inner) 'observation ', 'rule named ', method, call. = FALSE)
  label = if (inner) 'observation_options' else 'quadrature_options'
  ans = if (method == 'gl') {
      gl_options(if (inner) control$observation_nodes else control$nodes,
        if (inner) control$observation_options %||% list() else control$quadrature_options %||% list(),
        default_batch = if (inner) 128L else 8192L, label)
    } else staged_options(control, inner = TRUE)
  # Outer GH with a new PLN observation mesh keeps its established outer rule
  # (above). Panel score budgets take the evaluation's score floors (set at the
  # staged entry points); PLN observation rows check unit-free local derivatives.
  if (!inner) ans$score_floors = control$.score_floors
  ans
}

# Floors of the free coordinates' panel score budgets, tau*(floor + |score|).
staged_score_floor = function(options, free = TRUE) {
  if (is.null(options$score_floors)) 1 else options$score_floors[free]
}

staged_mesh_required = function(message, theta, diagnostics = NULL, probes = NULL) {
  if (is.null(probes)) probes = list(list(theta = theta, diagnostics = diagnostics))
  stop(staged_refinement_condition(paste('Staged quadrature:', message), diagnostics, probes))
}

# 2. Saved mesh bindings and compatibility -------------------------------

# A mesh is bound to its model, data, controls and restrictions. The data
# arrays enter through an MD5 digest of their raw little-endian bytes, storage
# types, dimensions and dimnames (no serialize() header, so the digest does not
# depend on the R version).
staged_binding = function(prepared, control, fixed) {
  data = list(x = prepared$x, w = prepared$w, y = prepared$y, d = prepared$d,
    active = prepared$active, offsets = list(prepared$offset_outcome, prepared$offset_selection),
    group = prepared$group)
  c(list(model = prepared$model_id %||% resolve_model(prepared$family,
      prepared$has_selection, prepared$structure)$id,
    structure = prepared$structure, family = prepared$family,
    outcome_scale = prepared$outcome_scale, outcome_prefix = prepared$outcome_prefix,
    has_selection = prepared$has_selection),
    list(data_md5 = staged_data_md5(data)),
    list(panels = prepared$n_panels, rows = prepared$n_rows,
    fixed = fixed, exact_reductions = control$exact_reductions, outer = control$quadrature,
    nodes = control$nodes, inner = control$observation_quadrature %||% 'gh',
    physical_pln = control$.physical_pln, observation_nodes = control$observation_nodes,
    quadrature_options = control$quadrature_options,
    observation_options = control$observation_options))
}

staged_data_md5 = function(data) {
  arrays = c(data[setdiff(names(data), 'offsets')],
    list(offset_outcome = data$offsets[[1L]], offset_selection = data$offsets[[2L]]))
  bytes = lapply(names(arrays), function(name) {
    value = arrays[[name]]
    header = paste(name, typeof(value), paste(dim(value) %||% length(value), collapse = 'x'),
      paste(c(unlist(dimnames(value)), names(value)), collapse = '\t'), sep = '\n')
    c(charToRaw(header), if (length(value)) writeBin(as.vector(value), raw(), endian = 'little'))
  })
  unname(tools::md5sum(bytes = unlist(bytes)))
}

staged_binding_matches = function(binding, prepared, control, fixed)
  identical(current_binding(binding), staged_binding(prepared, control, fixed))

# Legacy bindings hold the earlier control names (reduction, rule
# values adaptive, fixed and split, inner_nodes and one option list per rule
# with the earlier option names) and parameter names (the outcome scale and
# the restrictions); they are rewritten in the current layout.
current_binding = function(binding) {
  earlier = c('reduction', 'inner_nodes', 'split', 'gk', 'inner_split', 'inner_gk')
  # Bindings with current control fields may still use earlier parameter names.
  # Translate those names; the data digest and every data-derived label stay exact.
  if (any(names(binding$fixed) %in% names(renamed_parameters)))
    names(binding$fixed) = current_names(names(binding$fixed))
  if (is.character(binding$outcome_scale) && binding$outcome_scale %in% names(renamed_parameters))
    binding$outcome_scale = current_names(binding$outcome_scale)
  if (!any(earlier %in% names(binding))) return(binding)
  control = current_control(list(quadrature = binding$outer, inner_quadrature = binding$inner,
    split = binding$split, gk = binding$gk, inner_split = binding$inner_split,
    inner_gk = binding$inner_gk))
  head = binding[seq_len(match('fixed', names(binding)) - 1L)]
  if (is.character(head$outcome_scale)) head$outcome_scale = current_names(head$outcome_scale)
  fixed = binding$fixed
  if (!is.null(names(fixed))) names(fixed) = current_names(names(fixed))
  c(head, list(fixed = fixed, exact_reductions = binding$reduction,
    outer = control$quadrature, nodes = binding$nodes, inner = control$observation_quadrature,
    physical_pln = binding$physical_pln, observation_nodes = binding$inner_nodes,
    quadrature_options = control$quadrature_options,
    observation_options = control$observation_options))
}

# Legacy meshes may label GL cells 'split' and retain earlier parameter names
# and bindings. Normalize these labels before comparing mesh bindings;
# check_quadrature() reports whether refinement changed the mesh.
current_mesh = function(mesh) {
  # Current meshes occur on every evaluation. Parameter-only renames require
  # no walk through numerical cells; older rule layouts alone need relabelling.
  earlier = c('reduction', 'inner_nodes', 'split', 'gk', 'inner_split', 'inner_gk')
  old_rules = any(earlier %in% names(mesh$binding)) ||
    any(c(mesh$binding$outer, mesh$binding$inner) %in% names(old_rule_values))
  old_parameters = any(c(mesh$parameter_names, names(mesh$anchor),
    names(mesh$binding$fixed), mesh$binding$outcome_scale) %in% names(renamed_parameters))
  if (!old_rules && !old_parameters) return(mesh)
  relabel = function(x) {
    if (identical(x[['method', exact = TRUE]], 'split')) x$method = 'gl'
    for (i in seq_along(x)) if (is.list(x[[i]])) x[[i]] = relabel(x[[i]])
    x
  }
  for (part in if (old_rules) intersect(c('outer', 'inner'), names(mesh)) else character())
    if (is.list(mesh[[part]])) mesh[[part]] = relabel(mesh[[part]])
  if (is.character(mesh$parameter_names)) mesh$parameter_names = current_names(mesh$parameter_names)
  if (!is.null(names(mesh$anchor))) names(mesh$anchor) = current_names(names(mesh$anchor))
  if (is.list(mesh$binding)) mesh$binding = current_binding(mesh$binding)
  mesh
}

# 3. Positive cell rules and resource budgets ----------------------------

# A unit-cell tensor template depends only on the one-dimensional rule and
# the dimension. Templates are memoized here rather than stored with every
# panel or row mesh; meshes saved with their own template still use it.
staged_tensor_template = local({
  templates = list()
  function(rule, q) {
    if (!is.null(rule$tensor_template) && ncol(rule$tensor_template$nodes) == q)
      return(rule$tensor_template)
    # Rules name themselves (gl<order>); others use their exact values.
    key = paste(q, rule$key %||% paste(sprintf('%a', c(rule$nodes, rule$log_weights)),
      collapse = ' '))
    if (!is.null(templates[[key]])) return(templates[[key]])
    index = as.matrix(expand.grid(rep(list(seq_along(rule$nodes)), q)))
    ans = list(nodes = matrix(rule$nodes[index], nrow(index), q),
      log_weights = matrix(rule$log_weights[index], nrow(index), q))
    # Keep at most eight templates of 8192 numeric entries (64 KiB) each.
    # Higher supplementary orders use the same arithmetic without a copy.
    if (length(ans$nodes)*length(ans) <= 8192L) {
      if (length(templates) >= 8L) templates <<- tail(templates, 7L)
      templates[[key]] <<- ans
    }
    ans
  }
})

# The GL rule of order `nodes` per cell and its reference rule of order
# `reference_nodes` (above or below `nodes`).
staged_rule = function(nodes, reference_nodes, q = NULL) {
  rule = statmod::gauss.quad(nodes, 'legendre')
  reference = statmod::gauss.quad(reference_nodes, 'legendre')
  ans = list(nodes = rule$nodes, log_weights = log(rule$weights), reference =
    list(nodes = reference$nodes, log_weights = log(reference$weights),
      key = paste0('gl', reference_nodes)), key = paste0('gl', nodes))
  # A known dimension warms the shared template memo; meshes store only the
  # one-dimensional rules, each with its own template.
  if (!is.null(q) && q > 0L) {
    staged_tensor_template(ans, q)
    staged_tensor_template(ans[['reference']], q)
  }
  ans
}

staged_cells = function(axes) {
  q = length(axes)
  if (!q) return(list(lower = matrix(numeric(), 1L, 0L),
    upper = matrix(numeric(), 1L, 0L)))
  index = as.matrix(expand.grid(lapply(axes, function(x) seq_len(length(x)-1L))))
  lower = upper = matrix(0, nrow(index), q)
  for (j in seq_len(q)) {
    lower[, j] = axes[[j]][index[, j]]
    upper[, j] = axes[[j]][index[, j]+1L]
  }
  list(lower = lower, upper = upper)
}

staged_cell_rule = function(mesh, cells = seq_len(nrow(mesh$lower)), reference = FALSE) {
  q = ncol(mesh$lower); base = mesh$rule
  if (!q) return(list(nodes = matrix(numeric(), 1L, 0L), log_weights = 0,
    cells = 1L))
  if (reference) base = base[['reference']]
  template = staged_tensor_template(base, q)
  count = nrow(template$nodes); nc = length(cells)
  lower = mesh$lower[cells, , drop = FALSE]
  upper = mesh$upper[cells, , drop = FALSE]
  middle = (lower+upper)/2; half = (upper-lower)/2
  points = matrix(0, count*nc, q)
  for (j in seq_len(q)) points[, j] = rep(middle[, j], each = count)+
    rep(half[, j], each = count)*rep(template$nodes[, j], nc)
  list(nodes = points, log_weights = staged_cell_weights(mesh, cells, reference, template),
    cells = rep(cells, each = count))
}

staged_cell_weights = function(mesh, cells = seq_len(nrow(mesh$lower)), reference = FALSE,
                                template = NULL) {
  q = ncol(mesh$lower); base = mesh$rule
  if (!q) return(0)
  if (reference) base = base[['reference']]
  if (is.null(template)) template = staged_tensor_template(base, q)
  count = nrow(template$nodes); nc = length(cells)
  half = (mesh$upper[cells, , drop = FALSE]-mesh$lower[cells, , drop = FALSE])/2
  log_weights = numeric(count*nc)
  # Keep the original axis-wise addition order.
  for (j in seq_len(q))
    log_weights = log_weights+rep(log(half[, j]), each = count)+rep(template$log_weights[, j], nc)
  log_weights
}

staged_refine_cells = function(mesh, cells, dimension = NULL) {
  if (!length(cells)) return(mesh)
  if (is.null(dimension)) dimension = max.col(mesh$upper[cells, , drop = FALSE]-
    mesh$lower[cells, , drop = FALSE], ties.method = 'first')
  lower = mesh$lower[cells, , drop = FALSE]; upper = mesh$upper[cells, , drop = FALSE]
  midpoint = (lower+upper)/2
  left_upper = upper; right_lower = lower
  left_upper[cbind(seq_along(cells), dimension)] = midpoint[cbind(seq_along(cells), dimension)]
  right_lower[cbind(seq_along(cells), dimension)] = midpoint[cbind(seq_along(cells), dimension)]
  keep = setdiff(seq_len(nrow(mesh$lower)), cells)
  mesh$lower = rbind(mesh$lower[keep, , drop = FALSE], lower, right_lower)
  mesh$upper = rbind(mesh$upper[keep, , drop = FALSE], left_upper, upper)
  mesh
}

staged_extend_cells = function(mesh, bound) {
  # Extend a rectangle without replacing its already refined interior.
  q = ncol(mesh$lower); old = mesh$bound
  if (!q || bound <= old) return(mesh)
  extra = if (q == 1L) list(lower = matrix(c(-bound, old)),
    upper = matrix(c(-old, bound))) else list(
      lower = rbind(c(-bound, -bound), c(old, -bound), c(-old, -bound), c(-old, old)),
      upper = rbind(c(-old, bound), c(bound, bound), c(old, -old), c(old, bound)))
  mesh$lower = rbind(mesh$lower, extra$lower)
  mesh$upper = rbind(mesh$upper, extra$upper)
  mesh$bound = bound
  mesh
}

staged_node_count = function(mesh) {
  if (identical(mesh$route, 'physical_pln')) return(staged_physical_node_count(mesh))
  if (identical(mesh$route, 'independent'))
    return(staged_node_count(mesh$selection)+staged_node_count(mesh$outcome))
  q = ncol(mesh$lower)
  if (!q) return(0L)
  nrow(mesh$lower)*(length(mesh$rule$nodes)^q+length(mesh$rule[['reference']]$nodes)^q)
}

staged_interval_count = function(mesh) {
  if (is.null(mesh)) return(0L)
  if (identical(mesh$route, 'physical_pln'))
    return(as.integer(nrow(mesh$u$lower)*nrow(mesh$v$lower)))
  if (identical(mesh$route, 'independent'))
    return(staged_interval_count(mesh$selection)+staged_interval_count(mesh$outcome))
  if (!ncol(mesh$lower)) 0L else nrow(mesh$lower)
}

# `name` is the option list holding the limits; it is evaluated only on failure.
staged_resource_check = function(mesh, options, clock, name = NULL) {
  clock()
  intervals = staged_interval_count(mesh); nodes = staged_node_count(mesh)
  if (intervals > options$max_cells || nodes > options$max_nodes) {
    option = function(field) if (is.null(name)) field else paste0(name, '$', field)
    integration_error(paste0('Staged quadrature node or cell budget exhausted: ', intervals,
      ' cells (', option('max_cells'), ' = ', options$max_cells, '), ', format(nodes),
      ' nodes (', option('max_nodes'), ' = ', format(options$max_nodes), ')'))
  }
}

staged_inner_budget = function(prepared, options,
                                now = function() proc.time()[['elapsed']]) {
  # This is transient bookkeeping, never retained mesh state. Observation
  # kernels/refinement for all rows of one person share one cumulative limit;
  # outer-only work between calls is not charged to this inner budget.
  state = new.env(parent = emptyenv())
  state$panels = numeric(prepared$n_panels); state$total = 0
  state$active = FALSE; state$panel = NULL; state$row = NULL; state$started = 0
  snapshot = function(elapsed = NULL) {
    if (is.null(elapsed)) elapsed = if (state$active) max(0, now()-state$started) else 0
    panels = state$panels
    if (state$active) panels[state$panel] = panels[state$panel]+elapsed
    list(total_seconds = state$total+elapsed, panel_seconds = panels,
      seconds_total = options$seconds_total, seconds_per_panel = options$seconds_per_panel,
      active_row = if (state$active) state$row else NULL,
      scope = 'cumulative active inner wall time per person and top-level operation')
  }
  check = function() {
    elapsed = if (state$active) max(0, now()-state$started) else 0
    reason = if (state$total+elapsed > options$seconds_total) 'total' else
      if (!is.null(state$panel) && state$panels[state$panel]+elapsed > options$seconds_per_panel)
        'per-panel' else NULL
    if (!is.null(reason)) {
      option = if (reason == 'total') 'seconds_total' else 'seconds_per_panel'
      stop(structure(list(message = paste0('PLN observation quadrature ', reason,
        ' time budget exhausted (observation_options$', option, ' = ', options[[option]], ' s)'),
        call = NULL, diagnostics = c(snapshot(elapsed), list(panel = state$panel, row = state$row))),
        class = c('panelselect_inner_time_error', 'panelselect_integration_error', 'error', 'condition')))
    }
    invisible(NULL)
  }
  begin = function(row) {
    if (state$active) integration_error('PLN inner time accounting has an overlapping row')
    if (length(row) != 1L || !is.finite(row) || row < 1L || row > length(prepared$group))
      integration_error('PLN inner time accounting has no original row mapping')
    state$row = row; state$panel = prepared$group[row]
    state$started = now(); state$active = TRUE
    check()
    invisible(NULL)
  }
  end = function(check_limit = TRUE) {
    if (!state$active) return(invisible(NULL))
    elapsed = max(0, now()-state$started)
    state$panels[state$panel] = state$panels[state$panel]+elapsed
    state$total = state$total+elapsed; state$active = FALSE
    if (check_limit) check()
    invisible(NULL)
  }
  list(begin = begin, check = check, end = end, snapshot = snapshot)
}

staged_outer_batch = function(mesh, prepared, options, score_width) {
  q = ncol(mesh$lower); budget = staged_buffer_bytes
  active = prepared$active %||% rep(TRUE, prepared$n_rows)
  rows = max(c(0L, tabulate(prepared$group[active], prepared$n_panels)))
  if (!q) return(list(cells = 1L, requested_nodes = options$batch_nodes,
    effective_nodes = c(canonical = 0L, reference = 0L),
    estimated_buffer_bytes = 0, buffer_budget_bytes = budget,
    one_cell_over_budget = FALSE, worst_panel_rows = rows,
    score_width = score_width))
  canonical = length(mesh$rule$nodes)^q
  reference = length(mesh$rule[['reference']]$nodes)^q
  # A modest allocation estimate covers simultaneous high/reference node,
  # score and row-local buffers. Inner derivatives have their own 16 MiB cap.
  # Every cell retains its complete rule, even if a single cell exceeds this
  # estimate; the exception is explicit in diagnostics rather than truncating.
  entries = 8*(1+as.double(score_width)+q)+16*as.double(rows)
  bytes_per_cell = 8*(canonical+reference)*entries
  cells = min(nrow(mesh$lower), max(1L, floor(options$batch_nodes/
    max(canonical, reference))), max(1L, floor(budget/bytes_per_cell)))
  list(cells = as.integer(cells), requested_nodes = options$batch_nodes,
    effective_nodes = c(canonical = cells*canonical, reference = cells*reference),
    estimated_buffer_bytes = cells*bytes_per_cell, buffer_budget_bytes = budget,
    one_cell_over_budget = bytes_per_cell > budget, worst_panel_rows = rows,
    score_width = score_width)
}

staged_outer_radius = function(mesh) {
  if (identical(mesh$route, 'physical_pln')) return(max(
    staged_outer_radius(mesh$u), staged_outer_radius(mesh$v)))
  if (identical(mesh$route, 'independent')) return(max(
    staged_outer_radius(mesh$selection), staged_outer_radius(mesh$outcome)))
  if (!ncol(mesh$lower)) return(0)
  extrema = range(c(mesh$rule$nodes, mesh$rule[['reference']]$nodes))
  middle = (mesh$lower+mesh$upper)/2; half = (mesh$upper-mesh$lower)/2
  max(abs(middle+half*extrema[1L]), abs(middle+half*extrema[2L]))
}

# 4. Model layout, seed meshes and analytic tails ------------------------

staged_panel_layout = function(theta, prepared, control, fixed) {
  active = prepared$active %||% rep(TRUE, prepared$n_rows)
  if (prepared$has_selection && !any(active & prepared$d == 1))
    return(list(route = 'selection_only', parts = list(selection = selection_component(prepared))))
  independent = prepared$has_selection && !identical(prepared$structure, 'cross_sectional') &&
    prepared$family %in% c('probit','poisson','pln') && isTRUE(control$exact_reductions) &&
    is_fixed_zero('cor_re', theta, fixed) &&
    (prepared$family == 'poisson' || is_fixed_zero('cor_error', theta, fixed))
  physical = !independent && staged_physical_pln_eligible(theta, prepared, control, fixed)
  if (!independent && !physical) return(list(route = 'joint', parts = list(joint = prepared)))
  outcome = prepared
  outcome$has_selection = FALSE; outcome$model_id = NULL
  outcome$active = active & prepared$d == 1
  outcome$d = rep(1, prepared$n_rows); outcome$selected = outcome$active
  outcome$model_data = outcome$integration_cache = NULL
  list(route = if (physical) 'physical_pln' else 'independent',
    parts = list(selection = selection_component(prepared), outcome = outcome))
}

staged_seed_mesh = function(context, prepared, theta, control, fixed, inner = FALSE,
                             row = NULL) {
  options = staged_options(control, inner)
  method = if (inner) control$observation_quadrature else staged_outer_method(control)
  q = if (inner) 1L else context$q
  nodes = if (inner) control$observation_nodes else control$nodes
  rule = staged_rule(nodes, options$reference_nodes, q)
  axes = rep(list(seq(-options$bound, options$bound, by = if (q == 2L) 4 else 2)), q)
  hints = function(location, loading, target) {
    if (!q) return(invisible(NULL))
    for (j in seq_len(q)) {
      good = abs(loading[, j]) > sqrt(.Machine$double.eps)
      cuts = as.vector(outer(target, location[good], '-')/
        rep(loading[good, j], each = length(target)))
      axes[[j]] <<- sort(unique(c(axes[[j]], cuts[is.finite(cuts) &
        cuts > -options$bound & cuts < options$bound], options$bound)))
    }
    invisible(NULL)
  }
  if (inner) {
    B = context$b[row]; y = prepared$y[row]
    hints(B, matrix(context$error, 1L), log(pmax(y, .5)))
    if (prepared$has_selection && !context$tau_zero)
      hints(context$a[row], matrix(context$tau, 1L), 0)
  } else if (q) {
    if (prepared$has_selection) hints(context$a+context$delta*context$u$location[1L],
      matrix(rep(context$delta*context$u$loading[1L, ], each = prepared$n_rows),
        prepared$n_rows, q), if (q == 1L) c(-4,-2,-1,0,1,2,4) else 0)
    active = which((prepared$active %||% rep(TRUE, prepared$n_rows)) &
      (!prepared$has_selection | prepared$d == 1))
    if (length(active) && prepared$family != 'gaussian') {
      target = if (prepared$family == 'probit') if (q == 1L) c(-4,-2,-1,0,1,2,4) else 0 else
        log(pmax(prepared$y[active], .5))
      # Count peaks have row-specific targets; do not create all cross-pairs.
      location = context$b[active]+context$lambda*context$v$location[1L]
      loading = matrix(rep(context$lambda*context$v$loading[1L, ], each = length(active)),
        length(active), q)
      if (prepared$family == 'probit') hints(location, loading, target) else {
        for (k in seq_along(active)) hints(location[k], loading[k, , drop = FALSE], target[k])
      }
    }
  }
  cells = staged_cells(axes)
  c(cells, list(bound = options$bound, method = method, rule = rule))
}

staged_count_bounds = function(y) {
  maximum = exp(ifelse(y == 0, 0, y*log(pmax(y, 1))-y-lgamma(y+1)))
  first = y*maximum+exp((y+1)*log(y+1)-(y+1)-lgamma(y+1))
  list(mass = maximum, first = first)
}

# Log union-bound moments for omitted standard-normal coordinates. These bound
# the probability factor, not a Gaussian outcome density; derivative envelopes
# multiply the same moments before dividing by the calculated likelihood.
staged_tail_logs = function(bound, q) {
  if (!q) return(c(mass = -Inf, first = -Inf))
  mass = log(2)+pnorm(-bound, log.p = TRUE)
  first = log(2)+dnorm(bound, log = TRUE)
  if (q > 1L) first = logsum_weights(c(first,
    log(q-1)+.5*log(2/pi)+mass), 1L)$value
  c(mass = log(q)+mass, first = first)
}

staged_index_coefficients = function(theta, context, prepared, rows, kind) {
  q = context$q; p = length(theta); r = length(rows)
  map = if (kind == 'a') context$u else context$v
  scale = if (kind == 'a') context$delta else context$lambda
  scale_name = if (kind == 'a') 'sd_re_selection' else prepared$outcome_scale %||% 'sd_re_outcome'
  names_index = if (kind == 'a') context$data$selection_names else context$data$outcome_names
  design = if (kind == 'a') prepared$w[rows, , drop = FALSE] else prepared$x[rows, , drop = FALSE]
  ans = array(0, c(r, q+1L, p))
  if (!r) return(ans)
  for (k in seq_len(p)) {
    ans[, 1L, k] = scale*map$location_gradient[1L, k]
    if (q) for (j in seq_len(q)) ans[, j+1L, k] = scale*map$loading_gradient[1L, j, k]
    if (names(theta)[k] == scale_name) {
      ans[, 1L, k] = ans[, 1L, k]+map$location[1L]
      if (q) for (j in seq_len(q)) ans[, j+1L, k] = ans[, j+1L, k]+map$loading[1L, j]
    }
    index = match(names(theta)[k], names_index)
    if (!is.na(index)) ans[, 1L, k] = ans[, 1L, k]+design[, index]
  }
  ans
}

staged_outer_tail = function(theta, prepared, context, bound, loglik, score) {
  q = context$q; p = length(theta)
  logs = staged_tail_logs(bound, q)
  coefficients = matrix(0, p, q+1L)
  add = function(index, multiplier) {
    if (!length(index)) return(invisible(NULL))
    for (k in seq_len(p)) coefficients[k, ] <<- coefficients[k, ]+
      colSums(abs(matrix(index[, , k], nrow = dim(index)[1L]))*multiplier)
    invisible(NULL)
  }
  active = which(prepared$active %||% rep(TRUE, prepared$n_rows))
  selected = active[!prepared$has_selection | prepared$d[active] == 1]
  no = if (prepared$has_selection) active[prepared$d[active] == 0] else integer()
  phi0 = dnorm(0); c_tau = sqrt(1-context$tau^2)
  log_envelope = 0
  if (prepared$has_selection) add(staged_index_coefficients(theta, context, prepared, no, 'a'), phi0)
  B = staged_index_coefficients(theta, context, prepared, selected, 'b')
  if (prepared$family == 'gaussian') {
    A = staged_index_coefficients(theta, context, prepared, selected, 'a')
    sigma = context$error; tau = context$tau
    if (context$gaussian) {
      log_envelope = context$marginal$value[1L]
      coefficients[, 1L] = coefficients[, 1L]+abs(context$marginal$gradient[1L, ])
    } else {
      log_envelope = -length(selected)*log(sigma*sqrt(2*pi))
      add(B, exp(-.5)/sigma)
      k = match('sd_error', names(theta)); if (!is.na(k))
        coefficients[k, 1L] = coefficients[k, 1L]+length(selected)/sigma
    }
    if (length(selected)) {
      # Selected conditional indices are affine in the independent normals;
      # differentiate their coefficients before applying |Phi'| <= phi(0).
      av = cbind(context$a[selected]+context$delta*context$u$location[1L],
        matrix(rep(context$delta*context$u$loading[1L, ], each = length(selected)), length(selected), q))
      bv = cbind(context$b[selected]+context$lambda*context$v$location[1L],
        matrix(rep(context$lambda*context$v$loading[1L, ], each = length(selected)), length(selected), q))
      residual = -bv/sigma; residual[, 1L] = residual[, 1L]+prepared$y[selected]/sigma
      zeta = (av+tau*residual)/c_tau
      derivative = A/c_tau-tau*B/(sigma*c_tau)
      k = match('sd_error', names(theta)); if (!is.na(k))
        derivative[, , k] = derivative[, , k]-tau*residual/(sigma*c_tau)
      k = match('cor_error', names(theta)); if (!is.na(k))
        derivative[, , k] = derivative[, , k]+residual/c_tau+tau*zeta/c_tau^2
      add(derivative, phi0)
    }
  } else {
    A = if (prepared$has_selection) staged_index_coefficients(theta, context, prepared, selected, 'a') else NULL
    if (prepared$family == 'probit') {
      add(A, phi0); add(B, phi0)
      k = match('cor_error', names(theta)); if (!is.na(k))
        coefficients[k, 1L] = coefficients[k, 1L]+length(selected)/(2*pi*c_tau)
    } else {
      bounds = staged_count_bounds(prepared$y[selected])
      add(B, bounds$first)
      if (prepared$family == 'pln' && !is.null(context$inner)) {
        add(A, bounds$mass*phi0/c_tau)
        k = match('sd_error', names(theta)); if (!is.na(k))
          coefficients[k, 1L] = coefficients[k, 1L]+sum(bounds$first)*sqrt(2/pi)
        k = match('cor_error', names(theta)); if (!is.na(k))
          coefficients[k, 1L] = coefficients[k, 1L]+sum(bounds$mass)*
            (sqrt(2/pi)*phi0/c_tau+abs(context$tau)*dnorm(1)/c_tau^2)
      } else add(A, phi0)
    }
  }
  relative = exp(log_envelope+logs[['mass']]-loglik)
  score_tail = coefficients[, 1L]*relative+
    rowSums(coefficients[, -1L, drop = FALSE])*exp(log_envelope+logs[['first']]-loglik)+
    abs(score)*relative
  list(relative = relative, score = score_tail, log_envelope = log_envelope,
       derivative_coefficients = coefficients)
}

# 5. Cell reductions and evaluation records ------------------------------

# Node order that groups each cell's nodes for one batch reduction.
staged_cell_order = function(cells) {
  nc = length(unique(cells)); per_cell = length(cells)/nc
  list(cells = nc, reorder = as.vector(t(matrix(seq_along(cells), per_cell, nc))))
}

staged_reduce_cells = function(term, log_weights, cells, p, order = staged_cell_order(cells)) {
  reorder = order$reorder
  logsum_score(term$value[reorder]+log_weights[reorder],
    term$gradient[reorder, seq_len(p), drop = FALSE], order$cells)
}

staged_score_discrepancy = function(canonical_terms, reference_terms, weight, reference_weight,
                                    canonical_score, reference_score) {
  # Cellwise first-numerator differences expose cancellation, but normalizing
  # both rules by the canonical mass need not enclose their normalized scores.
  # Retain that diagnostic and also check the actual finite-rule score gap.
  numerator = colSums(abs(canonical_terms*weight-reference_terms*reference_weight))
  mass = sum(abs(weight-reference_weight))
  pmax(numerator+abs(canonical_score)*mass,
    abs(canonical_score-reference_score))
}

staged_recorder = function(mesh, theta) {
  ans = new.env(parent = emptyenv()); n = mesh$binding$rows
  ans$parameter_names = names(theta)
  ans$fixed_names = names(mesh$binding$fixed)
  ans$has_selection = mesh$binding$has_selection
  ans$selection_coordinate = as.integer(isTRUE(ans$has_selection) &&
    isTRUE(mesh$anchor['sd_re_selection'] > 0))
  ans$outcome_coordinate = if (isTRUE(mesh$anchor[
    mesh$binding$outcome_scale %||% 'sd_re_outcome'] > 0)) ans$selection_coordinate+1L else 0L
  ans$outer_radius = 0
  ans$outer_mass = numeric(mesh$binding$panels)
  ans$outer_moment = matrix(0, mesh$binding$panels, 2L)
  ans$factor_mass = ans$factor_moment = matrix(0, mesh$binding$panels, 2L)
  ans$factor_seen = matrix(FALSE, mesh$binding$panels, 2L)
  ans$inner_mass = rep(1, n)
  ans$inner_first = rep(sqrt(2/pi), n)
  ans$inner_work = numeric(n)
  ans$inner_error = rep(-Inf, n); ans$inner_tail = numeric(n)
  ans$inner_score_error = matrix(-Inf, n, length(theta))
  ans$inner_probes = vector('list', n)
  ans$inner_score_probes = vector('list', n)
  ans$inner_score_probe_error = rep(-Inf, n)
  ans
}

# 6. Outer evaluation and accuracy gates ---------------------------------

# Finite outer mass, absolute first moments and node radius of one batch.
staged_outer_measure = function(prepared, z, log_weights, shared = FALSE) {
  n = prepared$n_panels; q = ncol(z)
  weights = exp(log_weights)
  if (shared) {
    mass = rep(sum(weights), n)
    moment = vapply(seq_len(q), function(j) sum(weights*abs(z[, j])), numeric(1L))
    moment = matrix(rep(moment, each = n), n, q)
  } else {
    mass = rowSums(matrix(weights, n))
    moment = matrix(0, n, q)
    for (j in seq_len(q)) moment[, j] = rowSums(matrix(weights*abs(z[, j]), n))
  }
  list(radius = max(abs(z), 0), mass = mass, moment = moment, q = q)
}

staged_record_outer_measure = function(context, prepared, z, log_weights, shared = FALSE,
                                       measure = NULL) {
  recorder = context$staged_recorder %||% context$inner$recorder
  if (is.null(recorder) || isTRUE(context$staged_reference) ||
      isTRUE(context$inner[['reference']])) return(invisible(NULL))
  ids = prepared$.staged_panel_ids %||% seq_len(prepared$n_panels)
  if (is.null(measure)) measure = staged_outer_measure(prepared, z, log_weights, shared)
  q = measure$q; mass = measure$mass; moment = measure$moment
  recorder$outer_radius = max(recorder$outer_radius, measure$radius)
  component = context$staged_component %||% 'joint'
  if (component == 'joint') {
    recorder$outer_mass[ids] = recorder$outer_mass[ids]+mass
    if (q) recorder$outer_moment[ids, seq_len(q)] =
      recorder$outer_moment[ids, seq_len(q), drop = FALSE]+moment
  } else {
    k = match(component, c('selection','outcome'))
    recorder$factor_mass[ids, k] = recorder$factor_mass[ids, k]+mass
    if (q) recorder$factor_moment[ids, k] = recorder$factor_moment[ids, k]+moment[, 1L]
    recorder$factor_seen[ids, k] = TRUE
    measured = recorder$factor_mass[ids, , drop = FALSE]
    measured[!recorder$factor_seen[ids, , drop = FALSE]] = 1
    recorder$outer_mass[ids] = measured[, 1L]*measured[, 2L]
    if (recorder$selection_coordinate)
      recorder$outer_moment[ids, 1L] = recorder$factor_moment[ids, 1L]*measured[, 2L]
    if (recorder$outcome_coordinate)
      recorder$outer_moment[ids, recorder$outcome_coordinate] =
        recorder$factor_moment[ids, 2L]*measured[, 1L]
  }
  invisible(NULL)
}

staged_nested_error = function(theta, prepared, context, loglik, score, recorder) {
  p = length(theta)
  if (is.null(context$inner) || is.null(context$inner$quadrature) ||
      context$inner$quadrature != 'gl')
    return(list(probability = 0, score = numeric(p), probes = NULL))
  rows = if (prepared$has_selection) context$data$selected$rows else context$data$active$rows
  ids = (prepared$.staged_rows %||% seq_len(prepared$n_rows))[rows]
  if (!length(ids)) return(list(probability = 0, score = numeric(p), tail = 0, probes = list()))
  panel_id = (prepared$.staged_panel_ids %||% 1L)[1L]
  component = context$staged_component %||% 'joint'
  if (component == 'joint') {
    measure = recorder$outer_mass[panel_id]
    moments = recorder$outer_moment[panel_id, seq_len(context$q)]
  } else {
    k = match(component, c('selection','outcome'))
    measure = recorder$factor_mass[panel_id, k]
    moments = if (context$q) recorder$factor_moment[panel_id, k] else numeric()
  }
  if (!is.finite(measure) || measure <= 0)
    integration_error('PLN nested error accounting has no positive outer finite measure')
  maxima = vapply(prepared$y[rows], function(y) staged_count_bounds(y)$mass, numeric(1L))
  log_ceiling = log(maxima)+log(recorder$inner_mass[ids])
  log_companion = sum(log_ceiling)-log_ceiling
  # Positive high/reference/true count factors are bounded by M_y times their
  # finite normal mass. A telescoping product keeps all companion factors.
  probability = sum(exp(recorder$inner_error[ids]+log_companion-loglik))*measure
  first = colSums(exp(recorder$inner_score_error[ids, , drop = FALSE]+
    log_companion-loglik))*measure
  tail_logs = vapply(seq_along(ids), function(k) {
    bounds = staged_count_bounds(prepared$y[rows[k]])
    log(bounds$mass)+staged_tail_logs(context$inner$meshes[[ids[k]]]$bound, 1L)[['mass']]
  }, numeric(1L))
  tail = sum(exp(tail_logs+log_companion-loglik))*measure
  # Product perturbations use absolute row errors. This exposes repeated row
  # errors and prevents opposite inner/outer changes from cancelling in a gate.
  envelope = staged_outer_tail(theta, prepared, context, 1, loglik, score)$derivative_coefficients
  # These are the actual positive finite-rule moments. Adaptive GH weights
  # need not have unit prior mass or the moments of an exact normal density.
  inflation = max(c(1, recorder$inner_mass[ids], recorder$inner_first[ids]/sqrt(2/pi)))
  sup_first = inflation/exp(min(c(0, log_ceiling)))*(envelope[, 1L]+if (context$q)
    drop(envelope[, -1L, drop = FALSE] %*% moments)/measure else 0)
  # This is the direct natural-score perturbation at retained outer nodes.
  # Adaptive GH also differentiates its mode, transform and Jacobian; its
  # complete high/reference score comparison below remains empirical.
  list(probability = probability,
    score = first+probability*(sup_first+abs(score)),
    tail = tail, probes = staged_inner_probes(recorder, ids),
    log_companion = setNames(log_companion, as.character(ids)), measure = measure)
}

# One panel's canonical/reference discrepancy plus its propagated observation-
# rule error estimate and analytic omitted-tail bound.
staged_panel_error = function(theta, prepared, control, fixed, i, value, reference_value,
                              gradient, reference_gradient, recorder) {
  part = subset_panels(prepared, i); context = make_model_context(theta, part, control, fixed)
  nested = staged_nested_error(theta, part, context, value, gradient, recorder)
  list(nested = nested, dimension = context$q,
    change = abs(value-reference_value)+nested$probability,
    tail = max(nested$tail %||% 0, staged_inner_model_tail(part, context, value)),
    error = abs(gradient-reference_gradient)+nested$score)
}

staged_inner_model_tail = function(prepared, context, loglik) {
  if (is.null(context$inner) || is.null(context$inner$quadrature) ||
      context$inner$quadrature != 'gl') return(0)
  rows = if (prepared$has_selection) context$data$selected$rows else context$data$active$rows
  if (!length(rows)) return(0)
  ids = (prepared$.staged_rows %||% seq_len(prepared$n_rows))[rows]
  log_maxima = vapply(prepared$y[rows], function(y) log(staged_count_bounds(y)$mass), numeric(1L))
  tail_logs = vapply(ids, function(id) staged_tail_logs(
    context$inner$meshes[[id]]$bound, 1L)[['mass']], numeric(1L))
  # The outer normal density integrates to one, selection factors are <=1,
  # and every true count factor is <=M_y. No finite inside mass is substituted.
  exp(sum(log_maxima)+logsum_weights(tail_logs, 1L)$value-loglik)
}

# Parameter-independent work for one batch of cells: canonical and reference
# rules, their reduction orders and the outer-measure sums.
staged_cell_piece = function(mesh, cells, q, prepared) {
  rule = staged_cell_rule(mesh, cells)
  reference = staged_cell_rule(mesh, cells, reference = TRUE)
  piece = list(rule = rule, reference = reference, order = staged_cell_order(rule$cells),
    reference_order = staged_cell_order(reference$cells),
    measure_weights = rule$log_weights-q/2*log(2*pi)-rowSums(rule$nodes^2)/2)
  piece$measure = staged_outer_measure(prepared, rule$nodes, piece$measure_weights)
  piece
}

# Evaluator-local store of cell pieces, bounded at 64 MiB of node data. Without
# a workspace (direct evaluation or mesh preparation) nothing is kept.
staged_workspace_pieces = function(workspace, key, count) {
  stored = if (is.environment(workspace)) workspace[[key]]
  list(get = function(k) stored[[k]], put = function(k, piece) {
    if (!is.environment(workspace)) return(piece)
    bytes = 8*length(piece$rule$nodes)*3
    if ((workspace$bytes %||% 0)+bytes > staged_buffer_bytes) return(piece)
    workspace$bytes = (workspace$bytes %||% 0)+bytes
    if (is.null(workspace[[key]])) workspace[[key]] = vector('list', count)
    workspace[[key]][[k]] = piece
    piece
  })
}

staged_outer_panel = function(theta, prepared, control, fixed, mesh, recorder,
                              clock, score = TRUE, outer_radius = NULL) {
  if (identical(mesh$route, 'physical_pln'))
    return(staged_physical_pln_panel(theta, prepared, control, fixed, mesh, recorder, clock, score))
  options = staged_options(control)
  staged_resource_check(mesh, options, clock, staged_options_name(control))
  if (is.null(outer_radius)) outer_radius = staged_outer_radius(mesh)
  # Shared inner tau=0 kernels do not receive U separately. Their rho error
  # pullback uses one conservative actual-node radius for this entire panel,
  # including both orders/independent components. It is independent of batch
  # grouping and applies only to these fixed outer meshes, not adaptive GH.
  recorder$outer_radius = outer_radius
  if (identical(mesh$route, 'independent')) {
    parts = staged_panel_layout(theta, prepared, control, fixed)$parts
    if (!identical(names(parts), c('selection','outcome')))
      staged_mesh_required('Independent retained mesh no longer matches this panel', theta)
    selection = staged_outer_panel(theta, parts$selection, control, fixed,
      mesh$selection, recorder, clock, score, outer_radius)
    outcome = staged_outer_panel(theta, parts$outcome, control, fixed,
      mesh$outcome, recorder, clock, score, outer_radius)
    joint_loglik = selection$loglik+outcome$loglik
    joint_score = selection$score+outcome$score
    joint_control = control; joint_control$.staged_clock = clock
    joint_control$.staged_recorder = recorder
    joint_context = make_model_context(theta, prepared, joint_control, fixed)
    active_bounds = c(if (selection$dimension) mesh$selection$bound,
      if (outcome$dimension) mesh$outcome$bound)
    joint_tail = staged_outer_tail(theta, prepared, joint_context,
      if (length(active_bounds)) min(active_bounds) else options$bound,
      joint_loglik, joint_score)
    inner_model_tail = staged_inner_model_tail(prepared, joint_context, joint_loglik)
    return(list(loglik = joint_loglik, score = joint_score,
      reference_score = selection$reference_score+outcome$reference_score,
      loglik_error = selection$loglik_error+outcome$loglik_error,
      score_error = selection$score_error+outcome$score_error+joint_tail$score,
      tail = joint_tail$relative+inner_model_tail,
      outer_tail = joint_tail$relative, outer_tail_score = joint_tail$score,
      inner_tail = inner_model_tail,
      nested_error = selection$nested_error+outcome$nested_error,
      nested_score_error = selection$nested_score_error+outcome$nested_score_error,
      inner_probes = c(selection$inner_probes, outcome$inner_probes),
      nodes = selection$nodes+outcome$nodes, intervals = selection$intervals+outcome$intervals,
      batch = list(selection = selection$batch, outcome = outcome$batch),
      outer_radius = outer_radius,
      dimension = selection$dimension+outcome$dimension,
      factor_dimensions = c(selection = selection$dimension, outcome = outcome$dimension),
      route = 'independent selection and outcome components', selection = selection, outcome = outcome))
  }
  if (identical(mesh$route, 'selection_only')) prepared = selection_component(prepared)
  control$.staged_clock = clock; control$.staged_recorder = recorder
  context = make_model_context(theta, prepared, control, fixed)
  q = context$q; p = if (score) length(theta) else 0L
  if (q != ncol(mesh$lower)) staged_mesh_required('Latent dimension differs from the retained mesh', theta)
  kernel = make_panel_kernel(theta, prepared, control, fixed,
    order = as.integer(score), parameters = score, context = context)
  # The reference rule has its own nodes; its kernel also evaluates any staged
  # observation rule with that rule's reference.
  c_reference = control; c_reference$.staged_reference = TRUE
  reference_kernel = NULL
  # Cell arrays are cells x (1 + p): column 1 holds each cell's log mass and
  # the others its normalized score (the score integral divided by its mass).
  nc = nrow(mesh$lower); canonical_cells = reference_cells = matrix(0, nc, 1L+p)
  batch_diagnostics = staged_outer_batch(mesh, prepared, options, p)
  batch = batch_diagnostics$cells
  starts = seq.int(1L, nc, by = batch)
  pieces = staged_workspace_pieces(prepared$staged_workspace,
    paste(prepared$.staged_panel_ids[1L], context$staged_component, batch, sep = ':'),
    length(starts))
  for (k in seq_along(starts)) {
    cells = seq.int(starts[k], min(nc, starts[k]+batch-1L))
    piece = pieces$get(k) %||% pieces$put(k, staged_cell_piece(mesh, cells, q, prepared))
    rule = piece$rule
    staged_record_outer_measure(context, prepared, rule$nodes, piece$measure_weights,
      measure = piece$measure)
    term = kernel$evaluate(rule$nodes, prior = TRUE)
    canonical_cells[cells, ] = staged_reduce_cells(term, rule$log_weights, rule$cells, p, piece$order)
    ref_rule = piece$reference
    if (is.null(reference_kernel)) reference_kernel = make_panel_kernel(theta, prepared,
      c_reference, fixed, order = as.integer(score), parameters = score)
    reference_term = reference_kernel$evaluate(ref_rule$nodes, prior = TRUE)
    reference_cells[cells, ] = staged_reduce_cells(reference_term, ref_rule$log_weights,
      ref_rule$cells, p, piece$reference_order)
    clock()
  }
  if (any(is.na(canonical_cells[, 1L])) || any(is.na(reference_cells[, 1L])))
    staged_mesh_required('Outer likelihood is not finite on the retained mesh', theta)
  canonical = logsum_score(canonical_cells[, 1L], canonical_cells[, -1L, drop = FALSE], 1L)[1L, ]
  reference = logsum_score(reference_cells[, 1L], reference_cells[, -1L, drop = FALSE], 1L)[1L, ]
  if (!is.finite(canonical[1L]) || !is.finite(reference[1L]) ||
      (score && any(!is.finite(c(canonical[-1L], reference[-1L])))))
    staged_mesh_required('Outer likelihood or analytical score is not finite', theta)
  weight = exp(canonical_cells[, 1L]-canonical[1L]); ref_weight = exp(reference_cells[, 1L]-canonical[1L])
  cell_error = abs(weight-ref_weight)
  score_error = if (score) staged_score_discrepancy(canonical_cells[, -1L, drop = FALSE],
    reference_cells[, -1L, drop = FALSE], weight, ref_weight, canonical[-1L],
    reference[-1L]) else numeric(length(theta))
  full_score = if (score) canonical[-1L] else numeric(length(theta))
  tail = staged_outer_tail(theta, prepared, context, mesh$bound, canonical[1L], full_score)
  nested = staged_nested_error(theta, prepared, context, canonical[1L], full_score,
    recorder)
  inner_tail = max(nested$tail %||% 0,
    staged_inner_model_tail(prepared, context, canonical[1L]))
  list(loglik = canonical[1L], score = full_score,
    reference_score = if (score) reference[-1L] else numeric(length(theta)),
    loglik_error = max(abs(canonical[1L]-reference[1L]), sum(cell_error))+nested$probability,
    score_error = score_error+tail$score+nested$score,
    tail = tail$relative+inner_tail,
    outer_tail = tail$relative, inner_tail = inner_tail,
    outer_tail_score = tail$score,
    cell_error = cell_error, nested_error = nested$probability, nested_score_error = nested$score,
    inner_companion_log = nested$log_companion, inner_outer_mass = nested$measure,
    inner_probes = nested$probes,
    nodes = staged_node_count(mesh), intervals = if (q) nc else 0L, dimension = q,
    batch = batch_diagnostics,
    outer_radius = outer_radius,
    route = if (identical(mesh$route, 'selection_only')) 'selection only' else 'joint panel kernel')
}

staged_accuracy_reasons = function(value, options, theta, fixed, score = TRUE,
                                    loglik_target = options$panel_loglik_tol,
                                    score_target = options$panel_score_tol) {
  reasons = character()
  if (!is.finite(value$tail) || value$tail > options$tail_tol)
    reasons = c(reasons, 'omitted-tail envelope unresolved')
  if (!is.finite(value$loglik_error) || value$loglik_error > loglik_target)
    reasons = c(reasons, 'likelihood refinement unresolved')
  free = !names(theta) %in% names(fixed)
  if (score && any(!is.finite(value$score_error[free]) |
      value$score_error[free] > score_target*(staged_score_floor(options, free)+abs(value$score[free]))))
    reasons = c(reasons, 'analytical score refinement unresolved')
  reasons
}

# Panel accumulators shared by the staged evaluation routes. Canonical and
# reference scores are natural-parameter gradients; error estimates combine
# their absolute discrepancies with nested observation-rule estimates/tails.
# Per-cell error vectors belong to a live evaluation; a saved fit keeps each
# panel's scalar summaries (including independent components).
staged_panel_summaries = function(diagnostics) {
  summarize = function(panel) {
    if (!is.list(panel)) return(panel)
    panel[intersect(names(panel), c('cell_error', 'cell_refinement_error', 'dimension_error'))] = NULL
    for (part in c('selection', 'outcome')) if (is.list(panel[[part]])) panel[[part]] = summarize(panel[[part]])
    panel
  }
  if (is.list(diagnostics$panels)) diagnostics$panels = lapply(diagnostics$panels, summarize)
  diagnostics
}

staged_panel_record = function(theta, n) {
  zero = matrix(0, n, length(theta), dimnames = list(NULL, names(theta)))
  list(values = numeric(n), gradients = zero, reference = zero, errors = zero,
    changes = numeric(n), tails = numeric(n), diagnostics = vector('list', n), native = NULL)
}

# Native GH outer rule, including its differentiated adaptive geometry and
# standalone/fixed-zero reductions. Only the observation quadrature is staged.
staged_evaluate_native = function(theta, prepared, control, fixed, score, recorder, clock) {
  out = staged_panel_record(theta, prepared$n_panels)
  native = control; native$.staged_evaluation = TRUE
  native$.staged_recorder = recorder; native$.staged_clock = clock
  # An evaluator keeps its canonical and reference native caches in its memo;
  # each evaluation re-attaches its own recorder, clock and budget.
  memo = prepared$integration_cache$memo
  native_prepared = function(control, key) {
    held = if (is.environment(memo)) memo[[key]]
    ready = prepare_integration(held %||% prepared, control)
    if (is.environment(memo)) memo[[key]] = ready
    ready
  }
  canonical = evaluate_likelihood(theta, native_prepared(native, 'staged_canonical'),
    native, fixed, score)
  native$.staged_reference = TRUE
  reference = evaluate_likelihood(theta, native_prepared(native, 'staged_reference'),
    native, fixed, score)
  out$values = canonical$loglik
  if (score) { out$gradients = canonical$score; out$reference = reference$score }
  for (i in seq_len(prepared$n_panels)) {
    panel = staged_panel_error(theta, prepared, control, fixed, i, out$values[i],
      reference$loglik[i], out$gradients[i, ], out$reference[i, ], recorder)
    nested = panel$nested
    out$changes[i] = panel$change; out$tails[i] = panel$tail; out$errors[i, ] = panel$error
    out$diagnostics[[i]] = list(dimension = panel$dimension, nested_error = nested$probability,
      nested_score_error = nested$score,
      inner_companion_log = nested$log_companion, inner_outer_mass = nested$measure,
      inner_probes = nested$probes, tail = out$tails[i], loglik_error = out$changes[i],
      score_error = out$errors[i, ], score = out$gradients[i, ])
  }
  out$native = canonical$diagnostics
  out
}

# A batched direct kernel handles every cross-sectional observation and
# fixed-zero person-scale panel. PLN may still have its row-level integral.
staged_evaluate_direct = function(theta, prepared, control, fixed, score, recorder, clock) {
  n = prepared$n_panels; out = staged_panel_record(theta, n)
  direct_control = control; direct_control$.staged_recorder = recorder
  direct_control$.staged_clock = clock
  context = make_model_context(theta, prepared, direct_control, fixed)
  z = matrix(numeric(), n, 0L)
  staged_record_outer_measure(context, prepared, z, rep(0, n))
  kernel = make_panel_kernel(theta, prepared, direct_control, fixed,
    order = as.integer(score), parameters = score, context = context)
  canonical = kernel$evaluate(z, prior = TRUE)
  inner_changed = !is.null(context$inner) && !is.null(context$inner$quadrature) &&
    context$inner$quadrature == 'gl'
  if (inner_changed) {
    direct_control$.staged_reference = TRUE
    reference = make_panel_kernel(theta, prepared, direct_control, fixed,
      order = as.integer(score), parameters = score)$evaluate(z, prior = TRUE)
  } else reference = canonical
  out$values = canonical$value
  if (score) { out$gradients[, ] = canonical$gradient; out$reference[, ] = reference$gradient }
  for (i in seq_len(n)) {
    # Without staged observation rules a cross-sectional row needs no context.
    panel = if (inner_changed) staged_panel_error(theta, prepared, control, fixed, i,
      out$values[i], reference$value[i], out$gradients[i, ], out$reference[i, ], recorder) else {
      nested = list(probability = 0, score = numeric(length(theta)), probes = NULL)
      list(nested = nested, change = abs(out$values[i]-reference$value[i])+nested$probability,
        tail = 0, error = abs(out$gradients[i, ]-out$reference[i, ])+nested$score)
    }
    nested = panel$nested
    out$changes[i] = panel$change; out$tails[i] = panel$tail; out$errors[i, ] = panel$error
    out$diagnostics[[i]] = list(dimension = 0L, nodes = 0L, intervals = 0L,
      route = if (inner_changed) 'direct outer kernel with PLN row integration' else 'analytic direct kernel',
      nested_error = nested$probability, nested_score_error = nested$score,
      inner_probes = nested$probes, tail = out$tails[i], loglik_error = out$changes[i],
      score_error = out$errors[i, ], score = out$gradients[i, ])
  }
  out
}

# Retained positive outer meshes, one panel at a time under a per-panel budget.
staged_evaluate_mesh = function(theta, prepared, control, fixed, score, mesh, options,
                                recorder, clock) {
  out = staged_panel_record(theta, prepared$n_panels)
  for (i in seq_len(prepared$n_panels)) {
    tick = proc.time()[['elapsed']]
    panel_clock = staged_clock(options$seconds_per_panel, paste0(
      'Staged quadrature per-panel time budget exhausted at panel ', i, ' (',
      staged_options_name(control), '$seconds_per_panel = ', options$seconds_per_panel, ' s)'), clock)
    value = staged_outer_panel(theta, subset_panels(prepared, i), control, fixed,
      mesh$outer[[i]], recorder, panel_clock, score)
    out$values[i] = value$loglik; out$gradients[i, ] = value$score
    out$reference[i, ] = value$reference_score
    out$errors[i, ] = value$score_error; out$changes[i] = value$loglik_error
    out$tails[i] = value$tail
    value$elapsed_seconds = proc.time()[['elapsed']]-tick
    out$diagnostics[[i]] = value
  }
  out
}

evaluate_staged_integration = function(theta, prepared, control, fixed = numeric(), score = TRUE) {
  mesh = control$integration_mesh
  if (is.null(mesh)) staged_mesh_required('An immutable integration mesh must be prepared first', theta)
  # An evaluator's data, controls and mesh are fixed: check the binding once.
  workspace = prepared$staged_workspace
  if (!isTRUE(workspace$bound)) {
    if (!staged_binding_matches(mesh$binding, prepared, control, fixed) ||
        !identical(mesh$parameter_names, names(theta)))
      integration_error('Retained integration mesh is not bound to this model, data, controls and restrictions')
    if (is.environment(workspace)) workspace$bound = TRUE
  }
  control = prepare_rule_control(control)
  control$.score_floors = score_floors(prepared, theta, workspace)
  prepared$.staged_rows = seq_len(prepared$n_rows)
  prepared$.staged_panel_ids = seq_len(prepared$n_panels)
  prepared = prepare_model_data(prepared)
  options = staged_options(control)
  new_inner = staged_new_inner(prepared, control, theta)
  if (new_inner && is.null(control$.staged_inner_budget))
    control$.staged_inner_budget = staged_inner_budget(prepared, staged_options(control, TRUE))
  recorder = staged_recorder(mesh, theta)
  clock = staged_clock(staged_total_seconds(control, options), paste0(
    'Staged quadrature total time budget exhausted (', staged_options_name(control),
    '$seconds_total = ', staged_total_seconds(control, options), ' s)'), control$.staged_parent_clock)
  panels = if (control$quadrature %in% c('gh', 'agh') && !uses_gl_quadrature(control))
    staged_evaluate_native(theta, prepared, control, fixed, score, recorder, clock) else
    if (model_dimension(theta, prepared, control, fixed) == 0L)
      staged_evaluate_direct(theta, prepared, control, fixed, score, recorder, clock) else
      staged_evaluate_mesh(theta, prepared, control, fixed, score, mesh, options, recorder, clock)
  values = panels$values; gradients = panels$gradients; reference = panels$reference
  errors = panels$errors; changes = panels$changes; tails = panels$tails
  diagnostics = panels$diagnostics; native_diagnostics = panels$native
  clock()
  if (any(!is.finite(values)) || (score && any(!is.finite(c(gradients, reference)))))
    integration_error('The staged likelihood or its analytical score is not finite')
  if (length(fixed)) {
    gradients[, names(fixed)] = reference[, names(fixed)] = errors[, names(fixed)] = 0
  }
  # The rows are actual panels: repeated histories contribute repeatedly.
  # No normalization removes failures or treats repeated reference errors as noise.
  aggregate_loglik = sum(changes)
  info = parameter_layout(prepared); info$fixed = info$name %in% names(fixed)
  aggregate_score = colSums(errors)*abs(parameter_jacobian(theta, info))
  adaptive_inner = identical(control$quadrature, 'agh') &&
    model_dimension(theta, prepared, control, fixed) > 0L && prepared$family == 'pln' &&
    theta[['sd_error']] > 0 && (control$observation_quadrature %||% 'gh') == 'gl'
  outer_checked = !control$quadrature %in% c('gh','agh') &&
    any(vapply(diagnostics, function(x) (x$nodes %||% 0) > 0, logical(1L)))
  inner_checked = prepared$family == 'pln' && sum(recorder$inner_work) > 0 &&
    (control$observation_quadrature %||% 'gh') == 'gl'
  result = list(loglik = values, score = gradients, reference_score = reference,
    diagnostics = list(method = control$quadrature,
      observation_quadrature = control$observation_quadrature,
      dimension = model_dimension(theta, prepared, control, fixed),
      hybrid = prepared$family == 'pln' && theta[['sd_error']] > 0 &&
        ((control$quadrature %in% c('gh','agh') &&
          (control$observation_quadrature %||% 'gh') == 'gl') ||
         (control$quadrature == 'gl' &&
          (control$observation_quadrature %||% 'gh') %in% c('gh','agh'))),
      outer_order_checked = outer_checked, inner_order_checked = inner_checked,
      mesh_epoch = mesh$epoch, topology = 'fixed within stage', panels = diagnostics,
      aggregate_loglik_error = aggregate_loglik, aggregate_score_error = aggregate_score,
      maximum_tail_relative_bound = max(c(0, tails)), native = native_diagnostics,
      finite_prior_mass = if (any(recorder$outer_mass > 0)) recorder$outer_mass else NULL,
      finite_prior_absolute_moment = if (any(recorder$outer_mass > 0)) recorder$outer_moment else NULL,
      finite_prior_coordinate_system = vapply(diagnostics, function(x)
        x$coordinate_system %||% 'independent standard-normal person coordinates', character(1L)),
      inner_kernel_nodes = sum(recorder$inner_work),
      inner_time = if (new_inner) control$.staged_inner_budget$snapshot() else NULL,
      node_budget_scope = 'per retained outer or row-inner positive mesh; nested work is recorded separately',
      nested_score_status = if (adaptive_inner)
        'empirical full-native high/reference scores; direct affine envelope excludes adaptive-geometry sensitivity' else
        'empirical high/reference scores with analytic direct-kernel tail terms',
      order_checked = outer_checked || inner_checked, reference_status = if (adaptive_inner)
        'empirical refinement with analytic inner probability tails; adaptive-geometry score sensitivity is not enclosed' else
        'empirical refinement plus analytic omitted-tail envelopes'))
  if (!isTRUE(control$.staged_preparing)) {
    reasons = lapply(diagnostics, staged_accuracy_reasons, options = options,
      theta = theta, fixed = fixed, score = score)
    bad = which(lengths(reasons) > 0L)
    aggregate_bad = aggregate_loglik > options$total_loglik_tol ||
      (score && any(aggregate_score > options$total_score_tol*staged_score_floor(options)))
    if (length(bad) || aggregate_bad) {
      if (!length(bad)) bad = which(changes > 0 | rowSums(errors) > 0)
      probes = list(list(theta = theta, panels = bad, reasons = reasons,
        inner_probes = staged_inner_probes(recorder)))
      # One stable message keeps rejection tallies by reason; the failing
      # panels and their reasons travel in the probes (see staged_stage_summary()).
      staged_mesh_required('Likelihood, score or tail accuracy requires a new mesh', theta,
        result$diagnostics, probes)
    }
  }
  result
}

validate_staged_integration = function(theta, prepared, control, fixed = numeric(),
                                      mesh = NULL, information = FALSE) {
  if (!is.null(mesh)) control$integration_mesh = mesh
  result = evaluate_staged_integration(theta, prepared, control, fixed, score = TRUE)
  list(status = 'PASS', information = information, diagnostics = result$diagnostics)
}
