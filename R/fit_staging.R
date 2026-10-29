# A mesh is immutable during optimization and information differentiation.
# Rebuilding it requires a new complete optimization; all work shares three stages.
staged_information_tolerance = function(control, family) {
  options = list()
  if (uses_gl_quadrature(control)) options = c(options, list(control$quadrature_options))
  # Observation rules integrate only PLN likelihoods; their tolerances do not
  # tighten the information gate of other models.
  if (identical(family, 'pln') && identical(control$observation_quadrature, 'gl'))
    options = c(options, list(control$observation_options))
  min(vapply(options, function(x) x$information_tol, numeric(1)))
}

# The one mesh-refinement condition; staged_mesh_required() adds its prefix.
staged_refinement_condition = function(message, diagnostics = NULL, probes = NULL) {
  structure(list(message = message, call = NULL, probes = probes, diagnostics = diagnostics),
    class = c('panelselect_mesh_required', 'panelselect_integration_error', 'error', 'condition'))
}

merge_staged_rejections = function(records) {
  result = list(total = 0L, by_class = integer(), by_reason = integer(), last = NULL)
  for (record in records) {
    current = record$rejections
    if (is.null(current)) next
    result$total = result$total + current$total
    for (field in c('by_class', 'by_reason')) for (key in names(current[[field]])) {
      if (!key %in% names(result[[field]])) result[[field]][key] = 0L
      result[[field]][key] = result[[field]][key] + current[[field]][key]
    }
    if (!is.null(current$last)) result$last = current$last
  }
  result
}

# One clause per stage: its failing phase and the panels and reasons recorded.
# Accuracy-gate failures keep failing panel indices and per-panel reasons in
# their probes; other conditions keep them in their diagnostics.
staged_stage_summary = function(records) paste(vapply(records, function(record) {
  cause = record$error$diagnostics; probe = record$error$probes[[1L]]
  panels = probe$panels %||% cause$panel %||% cause$panels
  if (!is.numeric(panels)) panels = NULL
  reasons = probe$reasons %||% cause$reasons
  if (is.list(reasons) && length(panels) && length(reasons) >= max(panels)) reasons = reasons[panels]
  reasons = unique(unlist(reasons))
  paste0('stage ', record$epoch, ' (', record$phase,
    if (length(panels)) paste0(', panel', if (length(panels) > 1L) 's', ' ',
      paste(utils::head(panels, 5L), collapse = ', '), if (length(panels) > 5L) ', ...'),
    if (length(reasons)) paste0(': ', paste(utils::head(reasons, 3L), collapse = '; ')), ')')
}, character(1)), collapse = '; ')

fit_staged_model = function(theta, info, prepared, control, fixed) {
  # A refit/update starts a new bounded campaign; saved meshes serve evaluation.
  control$integration_mesh = control$.rules = NULL
  records = list(); previous = probes = options = NULL; optimize_time = 0
  information_computed = control$vcov != 'opg'
  for (epoch in seq_len(3L)) {
    started = proc.time()[['elapsed']]; evaluator = result = final = NULL
    start_theta = theta
    information_error = NULL; phase = 'prepare'; preparing = TRUE
    attempt = tryCatch({
      control$integration_mesh = prepare_staged_integration(theta, prepared, control,
        fixed, probes = probes, previous = previous, epoch = epoch)
      preparing = FALSE
      eta = parameter_to_optimizer(theta, info)
      evaluator = fit_evaluator(theta, info, prepared, control, fixed)
      initial = evaluator$evaluate(eta, trial = FALSE)
      # Keep the initial optimizer metric fixed across mesh restarts.
      if (is.null(options)) options = optimizer_options(control, prepared, info, initial$optimizer_score)
      phase = 'optimization'; tick = proc.time()[['elapsed']]
      result = tryCatch(optimize_model(eta, evaluator, control, options),
        finally = { optimize_time = optimize_time + proc.time()[['elapsed']] - tick })
      # The final evaluation applies the endpoint accuracy gate on this mesh;
      # repeating it would recompute the same deterministic integral.
      phase = 'endpoint_validation'
      final = evaluator$evaluate(result$par, trial = FALSE)
      theta = final$theta
      information = matrix(NA_real_, length(eta), length(eta), dimnames = list(names(eta), names(eta)))
      if (information_computed) {
        phase = 'information'
        information = information_matrix(result$par, evaluator,
          on_error = function(message) information_error <<- message)
      }
      TRUE
    }, error = identity)
    rejected = if (!is.null(evaluator)) evaluator$rejections() else NULL
    # A stage interrupted by a mesh request has no optimizer counts, but its
    # evaluator still counts the likelihood evaluations it spent.
    evaluations = if (is.function(evaluator$evaluations)) evaluator$evaluations() else 0L
    records[[epoch]] = list(epoch = epoch, status = if (inherits(attempt, 'error')) 'FAILED' else 'PASS',
      phase = phase, elapsed_seconds = proc.time()[['elapsed']] - started, evaluations = evaluations,
      start = start_theta,
      endpoint = if (!is.null(final)) final$theta else NULL,
      optimizer_counts = if (!is.null(result)) result$counts else NULL,
      optimizer_code = if (!is.null(result)) result$convergence else NULL,
      rejections = rejected, error = if (inherits(attempt, 'error'))
        list(class = class(attempt), message = conditionMessage(attempt),
          diagnostics = attempt[['diagnostics']], probes = attempt[['probes']]) else NULL)
    if (!inherits(attempt, 'error')) {
      final$diagnostics$staging = list(status = 'PASS', stages_used = epoch,
        maximum_stages = 3L, records = records,
        objective = if (identical(control$quadrature, 'agh'))
          paste('PLN observation meshes are immutable within each optimization/information stage;',
            'native outer adaptive GH geometry remains parameter dependent and fully differentiated; no tolerance relaxation.') else
          'One immutable positive finite rule per optimization/information stage; no tolerance relaxation.')
      return(list(control = control, evaluator = evaluator, result = result, final = final,
        options = options, optimize_time = optimize_time, rejections = merge_staged_rejections(records),
        evaluations = sum(vapply(records, function(record) record$evaluations, integer(1))),
        information = information, information_computed = information_computed,
        information_error = information_error))
    }
    # A retry would repeat a failed preparation with the same parameters, probes
    # and previous mesh under tighter targets, so preparation failures are final.
    if (!inherits(attempt, 'panelselect_mesh_required') || preparing || epoch == 3L) {
      failed = attempt
      failed$diagnostics = list(cause = attempt[['diagnostics']], stage_records = records,
        stages_used = epoch, maximum_stages = 3L, phase = phase,
        rejections = merge_staged_rejections(records))
      if (inherits(attempt, 'panelselect_mesh_required') && epoch == 3L) {
        failed$message = paste0('Integration mesh unresolved after three total stages: ',
          conditionMessage(attempt), '. By stage: ', staged_stage_summary(records),
          '. Details: the condition\'s $diagnostics$stage_records')
        class(failed) = unique(c('panelselect_mesh_exhausted', class(failed)))
      }
      stop(failed)
    }
    best = if (!is.null(evaluator)) evaluator$best_valid() else NULL
    if (is.null(final) && !is.null(best)) theta = best$theta
    probes = attempt[['probes']]
    if (is.null(probes)) probes = list(list(theta = theta, phase = phase))
    previous = control$integration_mesh
  }
}
