# Preparing one staged GL stage: seed the
# outer meshes, refine native, joint and independent panels and the observation
# meshes until the configured likelihood, score and tail gates pass, or report
# a bounded refinement failure. Bind the accepted mesh for every evaluation of
# that stage (evaluate_staged_integration); refinement is not part of the
# optimizer's within-stage objective or score.

# 1. Refine meshes between optimization stages ---------------------------

staged_select_cells = function(error) {
  positive = which(is.finite(error) & error > 0)
  if (!length(positive)) return(seq_len(min(8L, length(error))))
  order = positive[order(error[positive], decreasing = TRUE)]
  # Refine the cells accounting for most observed disagreement, bounded per
  # pass. Zero-error background cells retain their rule and their work.
  count = which(cumsum(error[order]) >= .7*sum(error[order]))[1L]
  head(order, min(64L, max(1L, count)))
}

# Seed or reuse each panel's outer meshes: one joint mesh, independent
# components, or physical PLN axes.
staged_seed_outer = function(mesh, theta, prepared, control, fixed, seed_control,
                             previous, options) {
  for (i in seq_len(prepared$n_panels)) {
    part = subset_panels(prepared, i)
    layout = staged_panel_layout(theta, part, control, fixed)
    children = list()
    for (component in names(layout$parts)) {
      piece = layout$parts[[component]]
      ctx = make_model_context(theta, piece, seed_control, fixed)
      old = previous$outer[[i]]
      if (layout$route %in% c('independent','physical_pln')) {
        child = if (layout$route == 'physical_pln')
          if (component == 'selection') 'u' else 'v' else component
        old = if (identical(old$route, layout$route)) old[[child]] else NULL
      }
      if (!is.null(old) && (identical(old$route, 'independent') ||
          !identical(ncol(old$lower), ctx$q))) old = NULL
      current = old %||% staged_seed_mesh(ctx, piece, theta, control, fixed)
      current$route = if (identical(layout$route, 'selection_only')) 'selection_only' else NULL
      current$method = staged_outer_method(control)
      current$rule = staged_rule(control$nodes, options$reference_nodes, ctx$q)
      if (current$bound < options$bound) current = staged_extend_cells(current, options$bound)
      children[[component]] = current
    }
    mesh$outer[[i]] = if (identical(layout$route, 'independent'))
      c(list(route = 'independent'), children) else if (identical(layout$route, 'physical_pln'))
        list(route = 'physical_pln', u = children$selection, v = children$outcome) else children[[1L]]
  }
  mesh
}

staged_refine_native = function(mesh, at, prepared, control, fixed, options,
                                selected_panels, target_log, target_score, clock, tightening) {
  # One native call checks every panel in a pass. Re-evaluating the complete
  # dataset separately for every panel would make preparation quadratic in N.
  for (iteration in seq_len(staged_refinement_passes)) {
    clock(); current_control = control; current_control$integration_mesh = mesh
    current_control$.staged_preparing = TRUE
    result = evaluate_staged_integration(at, prepared, current_control, fixed, TRUE)
    reasons = lapply(selected_panels, function(i) staged_accuracy_reasons(
      result$diagnostics$panels[[i]], options, at, fixed, TRUE,
      loglik_target = target_log, score_target = target_score))
    bad = selected_panels[lengths(reasons) > 0L]
    if (!length(bad)) break
    for (i in bad) {
      value = result$diagnostics$panels[[i]]
      if (!length(Filter(Negate(is.null), value$inner_probes)))
        staged_mesh_required('Native GH/PLN mesh cannot be refined from this probe', at,
          list(panel = i, reasons = reasons))
      mesh$inner = staged_refine_inner_probes(mesh$inner, value$inner_probes, value,
        result$loglik[i], min(target_log, target_score), options$tail_tol,
        at, prepared, control, fixed, clock, tightening)
    }
    if (iteration == staged_refinement_passes) staged_mesh_required(paste0(
      'Native GH/PLN observation mesh refinement limit exhausted (', staged_refinement_passes,
      ' passes; panels ', paste(utils::head(bad, 5L), collapse = ', '), if (length(bad) > 5L) ', ...',
      ': ', paste(unique(unlist(reasons)), collapse = '; '), ')'), at,
      list(panels = bad, reasons = reasons))
  }
  mesh
}

# Refine one outer leaf (a joint mesh or one independent component) until its
# likelihood, score and tail gates pass.
staged_refine_leaf = function(mesh, at, prepared, control, fixed, i, component, pieces,
                              independent, options, local_options, local_log, local_score,
                              new_inner, panel_clock, tightening) {
  for (iteration in seq_len(staged_refinement_passes)) {
    panel_clock(); staged_resource_check(mesh$outer[[i]], options, panel_clock,
      staged_options_name(control))
    current_control = control; current_control$integration_mesh = mesh
    current_control$.staged_preparing = TRUE
    recorder = staged_recorder(mesh, at)
    leaf = if (independent) mesh$outer[[i]][[component]] else mesh$outer[[i]]
    value = staged_outer_panel(at, pieces[[component]], current_control, fixed, leaf,
      recorder, panel_clock, TRUE)
    reasons = staged_accuracy_reasons(value, local_options, at, fixed, TRUE,
      loglik_target = local_log, score_target = local_score)
    if (!length(reasons)) break
    free = !names(at) %in% names(fixed)
    score_budget = local_score*(staged_score_floor(options, free)+abs(value$score[free]))
    outer_score_tail = value$outer_tail_score %||% numeric(length(at))
    inner_score_tail = value$inner_tail_score %||% numeric(length(at))
    # Two individually admissible positive tail components can still exhaust
    # the same score budget together; their sum alone then makes interior
    # acceptance impossible.
    analytic_score_tail = outer_score_tail[free]+inner_score_tail[free]
    combined_score_tail_bad = !is.finite(analytic_score_tail) |
      analytic_score_tail > score_budget
    outer_score_tail_bad = any(!is.finite(outer_score_tail[free]) |
      outer_score_tail[free] > score_budget |
      (combined_score_tail_bad & outer_score_tail[free] > 0))
    inner_score_tail_bad = any(!is.finite(inner_score_tail[free]) |
      inner_score_tail[free] > score_budget |
      (combined_score_tail_bad & inner_score_tail[free] > 0))
    inner_score_bad = any(!is.finite(value$nested_score_error[free]) |
      value$nested_score_error[free] > score_budget/4)
    physical = identical(leaf$route, 'physical_pln')
    if (new_inner && inner_score_tail_bad) {
      # The physical route exposes the true-model inner score tail, which the
      # fixed observation domain bounds. General nested errors are handled
      # through their recorded probes below.
      piece = pieces[[component]]
      rows = which((piece$active %||% rep(TRUE, piece$n_rows)) &
        (!piece$has_selection | piece$d == 1))
      ids = piece$.staged_rows[rows]
      staged_bound_failure('PLN omitted-inner-score-tail envelope unresolved at maximum bound',
        staged_options(control, TRUE)$bound, staged_options_name(control, TRUE, bound = TRUE), at,
        list(panel = i, rows = ids, component = component,
          bound = vapply(ids, function(id) mesh$inner[[id]]$bound, numeric(1L)),
          score_tail = inner_score_tail, score_budget = score_budget,
          inner_probes = value$inner_probes))
    } else if (new_inner && (value$nested_error > local_log/4 ||
        value$inner_tail > local_options$tail_tol || inner_score_bad)) {
      mesh$inner = staged_refine_inner_probes(mesh$inner, value$inner_probes, value,
        value$loglik, min(local_log, local_score), local_options$tail_tol,
        at, prepared, control, fixed, panel_clock, tightening)
    } else if (value$tail > local_options$tail_tol || outer_score_tail_bad) {
      envelope = if (physical) {
        if (outer_score_tail_bad) 'Physical PLN omitted-outer-score-tail' else 'Physical PLN omitted-outer-tail'
      } else if (outer_score_tail_bad) 'Outer omitted-score-tail' else 'Outer omitted-tail'
      staged_bound_failure(paste(envelope, 'envelope unresolved at maximum bound'), options$bound,
        staged_options_name(control, bound = TRUE), at,
        list(panel = i, component = component,
          bound = if (physical) c(u = leaf$u$bound, v = leaf$v$bound) else leaf$bound,
          tail = value$tail, score_tail = outer_score_tail, score_budget = score_budget))
    } else if (physical) leaf = staged_physical_refine(leaf, value) else
      leaf = staged_refine_cells(leaf, staged_select_cells(value$cell_error))
    if (independent) mesh$outer[[i]][[component]] = leaf else mesh$outer[[i]] = leaf
    if (iteration == staged_refinement_passes) staged_mesh_required(paste0(
      'Outer mesh refinement limit exhausted (', staged_refinement_passes, ' refinements of panel ', i,
      if (!identical(component, 'joint')) paste0(' ', component), ': ', paste(reasons, collapse = '; '), ')'), at,
      list(panel = i, component = component, reasons = reasons))
  }
  mesh
}

# After both independent components pass their own gates, check the joint
# model's union/density envelope. Returns NULL when it passes, otherwise the
# mesh with refined observation probes for another component round.
staged_check_independent_joint = function(mesh, at, prepared, control, fixed, i, part,
                                          pieces, options, target_log, target_score,
                                          new_inner, panel_clock, tightening, factor_iteration) {
  # A sum of canonical component-relative tails is not a true joint-tail
  # enclosure. Check the original joint model's union/density envelope
  # after both finite component objectives have met their own gates.
  current_control = control; current_control$integration_mesh = mesh
  current_control$.staged_preparing = TRUE
  recorder = staged_recorder(mesh, at)
  combined = staged_outer_panel(at, part, current_control, fixed,
    mesh$outer[[i]], recorder, panel_clock, TRUE)
  reasons = staged_accuracy_reasons(combined, options, at, fixed, TRUE,
    loglik_target = target_log, score_target = target_score)
  if (!length(reasons)) return(NULL)
  free = !names(at) %in% names(fixed)
  inner_score_tail = combined$inner_tail_score %||% numeric(length(at))
  score_budget = target_score*(staged_score_floor(options, free)+abs(combined$score[free]))
  inner_score_tail_bad = any(!is.finite(inner_score_tail[free]) |
    inner_score_tail[free] > score_budget)
  nested_score_bad = any(!is.finite(combined$nested_score_error[free]) |
    combined$nested_score_error[free] > score_budget)
  changed = FALSE
  if (new_inner && (combined$inner_tail > options$tail_tol/2 || inner_score_tail_bad)) {
    rows = which((part$active %||% rep(TRUE, part$n_rows)) & part$d == 1)
    ids = part$.staged_rows[rows]
    staged_bound_failure(if (inner_score_tail_bad)
        'Joint PLN omitted-inner-score-tail envelope unresolved at maximum bound' else
        'Joint PLN omitted-inner-tail envelope unresolved at maximum bound',
      staged_options(control, TRUE)$bound, staged_options_name(control, TRUE, bound = TRUE), at,
      list(panel = i, rows = ids, tail = combined$inner_tail,
        bound = vapply(ids, function(id) mesh$inner[[id]]$bound, numeric(1L)),
        score_tail = inner_score_tail, score_budget = score_budget))
  } else if (new_inner && nested_score_bad) {
    # Independent generic kernels retain a combined refinement/tail row
    # envelope, not a separate tail-only panel score. Tighten the recorded
    # outcome probes using the existing absolute derivative budget; inner
    # preparation stops if its analytic tail is the obstruction.
    current = mesh$inner
    mesh$inner = staged_refine_inner_probes(mesh$inner, combined$inner_probes,
      combined$outcome, combined$outcome$loglik, min(target_log, target_score),
      options$tail_tol, at, prepared, control, fixed, panel_clock, tightening)
    changed = !identical(current, mesh$inner)
  }
  if (combined$outer_tail > options$tail_tol/2 ||
      any(combined$outer_tail_score[free] > score_budget/4)) {
    active_components = names(pieces)[vapply(names(pieces), function(name)
      ncol(mesh$outer[[i]][[name]]$lower) > 0L, logical(1L))]
    if (length(active_components))
      staged_bound_failure('Joint omitted-outer-tail envelope unresolved at maximum bound',
        options$bound, staged_options_name(control, bound = TRUE), at,
        list(panel = i, components = active_components,
          bound = vapply(active_components, function(name) mesh$outer[[i]][[name]]$bound, numeric(1L)),
          tail = combined$outer_tail, score_tail = combined$outer_tail_score))
  }
  if (!changed || factor_iteration == staged_refinement_passes)
    staged_mesh_required('Independent component validation requires a new mesh', at,
      list(panel = i, reasons = reasons, value = combined),
      list(list(theta = at, panels = i, reasons = reasons, inner_probes = combined$inner_probes)))
  mesh
}

# 2. Prepare the next immutable stage ------------------------------------

prepare_staged_integration = function(theta, prepared, control, fixed = numeric(),
                                      probes = NULL, previous = NULL, epoch = 1L) {
  if (!is.numeric(epoch) || length(epoch) != 1L || epoch < 1 || epoch > 3)
    integration_error('Integration mesh epoch must be one of three total stages')
  control$.score_floors = score_floors(prepared, theta)
  prepared = prepare_model_data(prepared)
  prepared$.staged_rows = seq_len(prepared$n_rows)
  prepared$.staged_panel_ids = seq_len(prepared$n_panels)
  options = staged_options(control); started = proc.time()[['elapsed']]
  native_outer = control$quadrature %in% c('gh', 'agh') && !uses_gl_quadrature(control)
  clock = staged_clock(staged_total_seconds(control, options), paste0(
    'Staged mesh preparation total time budget exhausted (', staged_options_name(control),
    '$seconds_total = ', staged_total_seconds(control, options), ' s)'))
  control$.staged_parent_clock = clock
  tightening = 10^(epoch-1L)
  mesh = list(binding = staged_binding(prepared, control, fixed), parameter_names = names(theta),
    anchor = theta, epoch = as.integer(epoch), outer = vector('list', prepared$n_panels),
    inner = vector('list', prepared$n_rows), preparation = NULL)
  new_inner = staged_new_inner(prepared, control, theta)
  if (new_inner)
    control$.staged_inner_budget = staged_inner_budget(prepared, staged_options(control, TRUE))
  seed_control = control; seed_control$.staged_seed = TRUE
  context = make_model_context(theta, prepared, seed_control, fixed)
  if (new_inner) mesh = staged_seed_inner(mesh, theta, prepared, control, fixed,
    context, previous, clock, tightening)
  if (!native_outer && (context$q > 0L || new_inner))
    mesh = staged_seed_outer(mesh, theta, prepared, control, fixed, seed_control,
      previous, options)
  if (!native_outer && !new_inner && context$q == 0L) {
    clock()
    mesh$preparation = list(elapsed_seconds = proc.time()[['elapsed']]-started,
      probes_checked = 1L, factor = tightening, outer_intervals = rep(0L, prepared$n_panels),
      inner_intervals = rep(0L, prepared$n_rows), route = 'analytic direct kernel')
    return(mesh)
  }
  if (!is.null(probes$theta)) probes = list(probes)
  points = c(list(list(theta = theta)), Filter(function(x) is.list(x) &&
    is.numeric(x$theta) && identical(names(x$theta), names(theta)), probes %||% list()))
  # No cache or row mesh is mutated during evaluation. The following local
  # copies are refinements between stages, including all recorded failed probes.
  for (point in points) {
    at = point$theta
    if (new_inner && !is.null(point$inner_probes)) for (probe in point$inner_probes) {
      if (!is.null(probe$row) && !is.null(mesh$inner[[probe$row]]))
        mesh$inner[[probe$row]] = staged_prepare_inner(at, prepared, control, fixed,
          mesh$inner[[probe$row]], probe, clock, tightening)
    }
    selected_panels = point$panels %||% seq_len(prepared$n_panels)
    target_log = min(options$panel_loglik_tol, options$total_loglik_tol/prepared$n_panels)/(2*tightening)
    maximum_scale = max(c(1, abs(at[intersect(c('sd_re_selection','sd_re_outcome','sd_error'), names(at))])))
    target_score = min(options$panel_score_tol, options$total_score_tol/
      (prepared$n_panels*maximum_scale))/(2*tightening)
    if (native_outer) {
      mesh = staged_refine_native(mesh, at, prepared, control, fixed, options,
        selected_panels, target_log, target_score, clock, tightening)
      next
    }
    for (i in selected_panels) {
      panel_clock = staged_clock(options$seconds_per_panel, paste0(
        'Staged mesh preparation per-panel time budget exhausted at panel ', i, ' (',
        staged_options_name(control), '$seconds_per_panel = ', options$seconds_per_panel, ' s)'), clock)
      part = subset_panels(prepared, i)
      independent = identical(mesh$outer[[i]]$route, 'independent')
      pieces = if (independent) staged_panel_layout(at, part, control, fixed)$parts else list(joint = part)
      local_options = options
      if (independent) local_options$tail_tol = options$tail_tol/4
      local_log = target_log/length(pieces); local_score = target_score/length(pieces)
      for (factor_iteration in seq_len(if (independent) staged_refinement_passes else 1L)) {
        for (component in names(pieces)) mesh = staged_refine_leaf(mesh, at, prepared,
          control, fixed, i, component, pieces, independent, options, local_options,
          local_log, local_score, new_inner, panel_clock, tightening)
        if (!independent) break
        checked = staged_check_independent_joint(mesh, at, prepared, control, fixed, i,
          part, pieces, options, target_log, target_score, new_inner, panel_clock, tightening,
          factor_iteration)
        if (is.null(checked)) break
        mesh = checked
      }
    }
  }
  clock()
  # Without latent person effects the direct route evaluates every panel; its
  # zero-dimension outer meshes only carried the preparation loop.
  if (context$q == 0L) mesh$outer = vector('list', prepared$n_panels)
  mesh$preparation = list(elapsed_seconds = proc.time()[['elapsed']]-started,
    probes_checked = length(points), factor = tightening,
    outer_intervals = vapply(mesh$outer, staged_interval_count, integer(1L)),
    inner_intervals = vapply(mesh$inner, staged_interval_count, integer(1L)),
    inner_time = if (new_inner) control$.staged_inner_budget$snapshot() else NULL)
  mesh
}
