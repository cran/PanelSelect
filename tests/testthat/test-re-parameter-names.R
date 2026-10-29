# Stored objects can have the current control layout but use sd_id_* and cor_id.
# Only parameter labels are rewritten here; data IDs, variables, calls, rules
# and every numerical value remain exactly those of the original object.
re_saved_names = function(x) {
  earlier = c(sd_re_selection = 'sd_id_selection', sd_re_outcome = 'sd_id_outcome',
    cor_re = 'cor_id')
  hit = x %in% names(earlier)
  x[hit] = earlier[x[hit]]
  x
}

re_saved_numbers = function(x, parameter_rows = TRUE) {
  if (is.list(x) && !is.data.frame(x)) {
    for (i in seq_along(x)) if (!is.null(x[[i]]))
      x[[i]] = re_saved_numbers(x[[i]], parameter_rows)
  }
  if (is.numeric(x)) {
    if (!is.null(names(x))) names(x) = re_saved_names(names(x))
    if (!is.null(dimnames(x))) for (axis in seq_along(dimnames(x)))
      if ((parameter_rows || axis > 1L) && !is.null(dimnames(x)[[axis]]))
        dimnames(x)[[axis]] = re_saved_names(dimnames(x)[[axis]])
    accuracy = attr(x, 'checked_accuracy', exact = TRUE)
    if (is.list(accuracy)) attr(x, 'checked_accuracy') = re_saved_numbers(accuracy)
  }
  x
}

re_saved_control = function(control) {
  if (!is.null(names(control$parscale)))
    names(control$parscale) = re_saved_names(names(control$parscale))
  mesh = control$integration_mesh
  if (is.list(mesh)) {
    mesh$parameter_names = re_saved_names(mesh$parameter_names)
    names(mesh$anchor) = re_saved_names(names(mesh$anchor))
    names(mesh$binding$fixed) = re_saved_names(names(mesh$binding$fixed))
    if (is.character(mesh$binding$outcome_scale))
      mesh$binding$outcome_scale = re_saved_names(mesh$binding$outcome_scale)
    control$integration_mesh = mesh
  }
  control
}

re_saved_fit = function(fit) {
  for (field in c('coefficients', 'fixed', 'covariance', 'optimizer'))
    fit[[field]] = re_saved_numbers(fit[[field]])
  fit$integration = re_saved_numbers(fit$integration, parameter_rows = FALSE)
  colnames(fit$score) = re_saved_names(colnames(fit$score))
  fit$parameter_info$name = re_saved_names(fit$parameter_info$name)
  if (is.character(fit$prepared$outcome_scale))
    fit$prepared$outcome_scale = re_saved_names(fit$prepared$outcome_scale)
  fit$control = re_saved_control(fit$control)
  fit
}

re_name_data = function() {
  data.frame(id = rep(c('cor_id', 'sd_id_outcome', 'sd_id_selection', 'panel4', 'panel5', 'panel6'), each = 2),
    sd_id_outcome = rep(c(-1, 1), 6), d = c(1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1),
    y = c(.7, 1.5, 1, NA, NA, 1.4, .4, .9, NA, NA, 1.1, 1.8),
    count = c(0, 2, 1, 3, 2, 0, 1, 1, 4, 2, 0, 1))
}

re_name_fit = function(rule = 'gh') {
  data = re_name_data()
  if (rule == 'gh') {
    fixed = c(sd_re_selection = .6, sd_re_outcome = .4, sd_error = .7,
      cor_re = .2, cor_error = .1)
    panel_select(stats::as.formula('y ~ 1', env = baseenv()),
      stats::as.formula('d ~ 1', env = baseenv()), data, 'id', 'gaussian', fixed = fixed,
      start = c('selection.(Intercept)' = .2, 'outcome.(Intercept)' = 1), init = 'pooled',
      vcov = 'observed', control = panel_control(nodes = 6L, maxit = 0L,
        stationarity_warning = NULL))
  } else {
    # A genuine retained GL mesh with a free random-effect SD and modern fields.
    suppressWarnings(panel_select(stats::as.formula('count ~ sd_id_outcome', env = baseenv()),
      NULL, data, 'id', 'poisson',
      start = c('outcome.(Intercept)' = .2, outcome.sd_id_outcome = .1, sd_re_outcome = .5),
      init = 'pooled', control = panel_control(quadrature = 'gl', nodes = 12L, maxit = 0L,
        stationarity_warning = NULL, parscale = c('outcome.(Intercept)' = 1,
          outcome.sd_id_outcome = 2, sd_re_outcome = 3))))
  }
}

test_that('stored parameter-only controls and mesh bindings normalize without changing rules', {
  fit = re_name_fit('gl')
  old = re_saved_control(fit$control)
  expect_identical(old$quadrature, 'gl')
  expect_identical(names(old), names(fit$control))
  expect_identical(old$integration_mesh$outer, fit$control$integration_mesh$outer)
  expect_true('sd_id_outcome' %in% names(old$parscale))
  expect_identical(current_control(old), fit$control)
  expect_identical(validate_panel_control(old), validate_panel_control(fit$control))
  expect_identical(current_mesh(old$integration_mesh), fit$control$integration_mesh)
  expect_identical(current_mesh(fit$control$integration_mesh), fit$control$integration_mesh)
  # A current-layout binding can itself contain old fixed/scale labels.
  binding = fit$control$integration_mesh$binding
  binding$fixed = c(sd_re_outcome = 0, cor_re = 0)
  binding$outcome_scale = 'sd_re_outcome'
  old_binding = binding
  names(old_binding$fixed) = re_saved_names(names(binding$fixed))
  old_binding$outcome_scale = 'sd_id_outcome'
  expect_identical(current_binding(old_binding), binding)
  expect_identical(old_binding$data_md5, binding$data_md5)
  expect_identical(current_names(c('sd_selection', 'sd_outcome', 'cor_id')),
    c('sd_re_selection', 'sd_re_outcome', 'cor_re'))
})

test_that('saved information metadata changes labels without changing numbers or caches', {
  fit = re_name_fit()
  information = fit$optimizer$information
  # Retained GL accuracy records use parameter-labelled matrices as attributes.
  # This independent metadata fixture isolates translation from differentiation.
  attr(information, 'checked_accuracy') = list(status = 'PASS',
    canonical_information = information, predecessor_information = information,
    native_optimizer_step = setNames(c(.001, .001), rownames(information)))
  fit$optimizer$information = information
  fit$integration$score_by_unit = fit$score
  cache = new.env(parent = emptyenv())
  cache$parameter_names = c('sd_id_selection', 'cor_id')
  fit$prepared$integration_cache = cache
  old = re_saved_fit(fit)
  # Include a structural parameter in an information check even when this
  # deliberately small fit fixes it; metadata must be translated at any depth.
  dimnames(old$optimizer$information) = rep(list(c('sd_id_selection', 'cor_id')), 2)
  dimnames(fit$optimizer$information) = rep(list(c('sd_re_selection', 'cor_re')), 2)
  old_accuracy = attr(old$optimizer$information, 'checked_accuracy')
  new_accuracy = attr(fit$optimizer$information, 'checked_accuracy')
  for (field in c('canonical_information', 'predecessor_information')) {
    dimnames(old_accuracy[[field]]) = dimnames(old$optimizer$information)
    dimnames(new_accuracy[[field]]) = dimnames(fit$optimizer$information)
  }
  attr(old$optimizer$information, 'checked_accuracy') = old_accuracy
  attr(fit$optimizer$information, 'checked_accuracy') = new_accuracy
  before = serialize(old, NULL)
  normalized = current_fit(old)
  expect_identical(normalized, fit)
  expect_identical(serialize(old, NULL), before)
  expect_identical(normalized$prepared$integration_cache, cache)
  expect_identical(cache$parameter_names, c('sd_id_selection', 'cor_id'))
  expect_identical(normalized$data, old$data)
  expect_identical(rownames(normalized$score), sort(unique(old$data$id)))
  expect_identical(rownames(normalized$integration$score_by_unit), rownames(old$score))
  expect_identical(names(normalized$loglik_by_unit), names(old$loglik_by_unit))
  expect_identical(current_fit(normalized), normalized)
})

test_that('stored parameter labels retain prediction, information, update and ID behavior', {
  for (rule in c('gh', 'gl')) {
    fit = re_name_fit(rule)
    old = re_saved_fit(fit)
    path = tempfile(fileext = '.rds')
    on.exit(unlink(path), add = TRUE)
    saveRDS(old, path, version = 2)
    before = tools::md5sum(path)
    old = readRDS(path)
    expect_identical(current_fit(old), fit)
    expect_identical(coef(old), coef(fit))
    expect_identical(vcov(old), vcov(fit))
    expect_identical(check_stationarity(old), check_stationarity(fit))
    expect_identical(utils::capture.output(print(old)), utils::capture.output(print(fit)))
    expect_identical(utils::capture.output(print(summary(old))),
      utils::capture.output(print(summary(fit))))
    expect_identical(predict(old, fit$data[1:4, ]), predict(fit, fit$data[1:4, ]))
    expect_identical(predict(old, fit$data[1:4, ], se.fit = TRUE),
      predict(fit, fit$data[1:4, ], se.fit = TRUE))
    expect_identical(update(old, evaluate = FALSE), update(fit, evaluate = FALSE))
    # Explicit init matters for genuinely older saved ProbitRE controls; this
    # parameter-name migration alone leaves supported stored settings intact.
    expect_identical(coef(suppressWarnings(update(old))), coef(suppressWarnings(update(fit))))
    if (rule == 'gh') {
      expect_identical(add_information(old), add_information(fit))
      expect_identical(random_effects(old), random_effects(fit))
      expect_identical(predict(old, fit$data[1:4, ], conditioning = 'history', history = fit$data),
        predict(fit, fit$data[1:4, ], conditioning = 'history', history = fit$data))
    }
    expect_identical(tools::md5sum(path), before)
    expect_identical(readRDS(path), old)
  }
})

test_that('printing a saved summary translates only structural parameter labels', {
  fit = re_name_fit()
  current = summary(fit)
  old = current
  rownames(old$coefficients) = re_saved_names(rownames(old$coefficients))
  old$fixed = re_saved_names(old$fixed)
  names(old$parameter_legend) = re_saved_names(names(old$parameter_legend))
  old$fixed_values = gsub('sd_re_', 'sd_id_', old$fixed_values, fixed = TRUE)
  old$fixed_values = gsub('cor_re', 'cor_id', old$fixed_values, fixed = TRUE)
  before = serialize(old, NULL)
  expect_identical(utils::capture.output(print(old)), utils::capture.output(print(current)))
  expect_identical(serialize(old, NULL), before)
})
