# These references use the retained positive finite nodes, but independently
# form the scalar selection and count products and their natural derivatives.
staged_reduction_fixture = function(family = 'poisson', panels = 1:3) {
  data = data.frame(id = rep(1:3, each = 3), d = c(0, 0, 0, 1, 0, 1, 1, 1, 0),
    y = c(NA, NA, NA, 0, NA, 2, 1, 3, NA),
    x = c(-.83, .27, .91, -.57, .41, 1.13, -.71, .19, .87),
    w = c(.23, -.47, .81, .37, -.59, .73, .29, -.43, .97),
    os = c(.043, -.031, .071, -.057, .023, .089, -.019, .037, -.061),
    oy = c(.021, -.017, .051, -.043, .033, .067, -.029, .047, -.013))
  if (family == 'probit') data$y[data$d == 1] = c(0, 1, 0, 1)
  if (family == 'gaussian') data$y[data$d == 1] = c(-.6, .2, .8, -.3)
  # A subset keeps the listed panels' rows exactly; others are not prepared.
  data = data[data$id %in% panels, , drop = FALSE]
  prepared = prepare_data(y ~ x + offset(oy), d ~ w + offset(os), data, 'id', family, 'fail')
  info = parameter_layout(prepared)
  theta = c('selection.(Intercept)' = .173, selection.w = -.283,
    'outcome.(Intercept)' = -.237, outcome.x = .193,
    sd_re_selection = .47, sd_re_outcome = .61, sd_error = .53, cor_re = -.35, cor_error = .25)[info$name]
  list(data = data, prepared = prepared, theta = theta)
}

staged_reduction_control = function() staged_test_budget(panel_control(quadrature = 'gl', nodes = 12L,
  quadrature_options = list(reference_nodes = 15L),
  observation_quadrature = 'gl', observation_nodes = 12L,
  observation_options = list(reference_nodes = 15L)))

staged_reduction_weights = function(mesh = NULL, normal = NULL) {
  if (!is.null(normal)) return(list(nodes = normal$nodes[, 1L], log_weights = normal$log_weights))
  rule = staged_cell_rule(mesh)
  if (!ncol(rule$nodes)) return(list(nodes = 0, log_weights = 0))
  stopifnot(ncol(rule$nodes) == 1L)
  list(nodes = rule$nodes[, 1L], log_weights = rule$log_weights + dnorm(rule$nodes[, 1L], log = TRUE))
}

staged_reduction_mix = function(logmass, score) {
  high = max(logmass); value = high + log(sum(exp(logmass - high)))
  list(value = value, score = colSums(score*exp(logmass - value)))
}

staged_selection_reference = function(theta, prepared, rows, weights) {
  names_alpha = paste0('selection.', colnames(prepared$w))
  a = drop(prepared$w[rows, , drop = FALSE] %*% theta[names_alpha]) + prepared$offset_selection[rows]
  signs = 2*prepared$d[rows] - 1
  z = sweep(sweep(outer(weights$nodes, rep(theta[['sd_re_selection']], length(rows))), 2L, a, '+'),
    2L, signs, '*')
  lp = pnorm(z, log.p = TRUE)
  partial = sweep(exp(dnorm(z, log = TRUE) - lp), 2L, signs, '*')
  score = matrix(0, length(weights$nodes), length(theta), dimnames = list(NULL, names(theta)))
  score[, names_alpha] = partial %*% prepared$w[rows, , drop = FALSE]
  score[, 'sd_re_selection'] = weights$nodes*rowSums(partial)
  staged_reduction_mix(weights$log_weights + rowSums(lp), score)
}

staged_count_reference = function(theta, prepared, rows, weights, inner_meshes = NULL) {
  names_beta = paste0('outcome.', colnames(prepared$x))
  b = drop(prepared$x[rows, , drop = FALSE] %*% theta[names_beta]) + prepared$offset_outcome[rows]
  value = numeric(length(weights$nodes))
  score = matrix(0, length(weights$nodes), length(theta), dimnames = list(NULL, names(theta)))
  for (j in seq_along(weights$nodes)) {
    v = weights$nodes[j]
    for (k in seq_along(rows)) {
      e_rule = if (prepared$family == 'pln') staged_reduction_weights(inner_meshes[[rows[k]]]) else
        list(nodes = 0, log_weights = 0)
      eta = b[k] + theta[['sd_re_outcome']]*v +
        if (prepared$family == 'pln') theta[['sd_error']]*e_rule$nodes else 0
      mu = exp(eta)
      row = staged_reduction_mix(dpois(prepared$y[rows[k]], mu, log = TRUE) + e_rule$log_weights,
        cbind(index = prepared$y[rows[k]] - mu,
          error = e_rule$nodes*(prepared$y[rows[k]] - mu)))
      value[j] = value[j] + row$value
      score[j, names_beta] = score[j, names_beta] + row$score[['index']]*prepared$x[rows[k], ]
      score[j, 'sd_re_outcome'] = score[j, 'sd_re_outcome'] + v*row$score[['index']]
      if (prepared$family == 'pln') score[j, 'sd_error'] = score[j, 'sd_error'] + row$score[['error']]
    }
  }
  staged_reduction_mix(weights$log_weights + value, score)
}

test_that('never-response staged panels retain the exact one-dimensional selection finite rule', {
  for (family in c('probit', 'gaussian', 'poisson', 'pln')) {
    # Only never-response panel 1 is checked; preparing panels 2-3 costs minutes.
    f = staged_reduction_fixture(family, panels = 1L); control = staged_reduction_control()
    control$integration_mesh = prepare_staged_integration(f$theta, f$prepared, control)
    mesh = control$integration_mesh$outer[[1L]]
    expect_identical(mesh$route, 'selection_only')
    expect_identical(ncol(mesh$lower), 1L)
    actual = evaluate_model(f$theta, f$prepared, control)
    reference = staged_selection_reference(f$theta, f$prepared,
      which(f$prepared$group == 1L), staged_reduction_weights(mesh))
    expect_equal(actual$loglik[1L], reference$value, tolerance = 2e-12)
    expect_equal(drop(actual$score[1L, ]), reference$score, tolerance = 2e-11)
    irrelevant = setdiff(names(f$theta), c(grep('^selection\\.', names(f$theta), value = TRUE), 'sd_re_selection'))
    expect_equal(actual$score[1L, irrelevant], rep(0, length(irrelevant)), ignore_attr = TRUE)
    expect_identical(actual$diagnostics$panels[[1L]]$dimension, 1L)
    expect_equal(actual$diagnostics$panels[[1L]]$nodes, staged_node_count(mesh))
    # Stored identity/layout remains that of the selected model, not its reduction.
    expect_identical(control$integration_mesh$parameter_names, names(f$theta))
    expect_identical(control$integration_mesh$binding$model, f$prepared$model_id)
  }
})

test_that('fixed zero person and observation correlations factor finite count rules with score scattering', {
  for (family in c('poisson', 'pln')) {
    f = staged_reduction_fixture(family)
    fixed_names = intersect(c('cor_re', 'cor_error'), names(f$theta))
    fixed = setNames(rep(0, length(fixed_names)), fixed_names)
    f$theta[names(fixed)] = fixed
    control = staged_reduction_control()
    control$integration_mesh = prepare_staged_integration(f$theta, f$prepared, control, fixed)
    actual = evaluate_model(f$theta, f$prepared, control, fixed)
    for (i in 2:3) {
      mesh = control$integration_mesh$outer[[i]]
      expect_identical(mesh$route, 'independent')
      rows = which(f$prepared$group == i); selected = rows[f$prepared$d[rows] == 1]
      s = staged_selection_reference(f$theta, f$prepared, rows, staged_reduction_weights(mesh$selection))
      y = staged_count_reference(f$theta, f$prepared, selected, staged_reduction_weights(mesh$outcome),
        control$integration_mesh$inner)
      reference = s$score + y$score; reference[names(fixed)] = 0
      expect_equal(actual$loglik[i], s$value + y$value, tolerance = 2e-12)
      expect_equal(drop(actual$score[i, ]), reference, tolerance = 2e-11)
      expect_equal(actual$diagnostics$panels[[i]]$factor_dimensions, c(selection = 1L, outcome = 1L))
      nodes = staged_node_count(mesh$selection) + staged_node_count(mesh$outcome)
      expect_equal(actual$diagnostics$panels[[i]]$nodes, nodes)
      expect_lt(nodes, staged_node_count(mesh$selection)*staged_node_count(mesh$outcome))
      # The joint true-model union tail cannot substitute finite component
      # probabilities for unknown exact inside masses. A ratio checks this
      # even when the absolute tail is far below ordinary test tolerances.
      bound = min(mesh$selection$bound, mesh$outcome$bound)
      joint_tail = 4*pnorm(-bound)*exp(-actual$loglik[i])
      expect_equal(actual$diagnostics$panels[[i]]$outer_tail/joint_tail,
        1, tolerance = 2e-12)
    }
    # A free rho at zero must retain the joint mesh and its surviving rho score.
    free = fixed[names(fixed) != 'cor_re']; joint = control; joint$integration_mesh = NULL
    joint$integration_mesh = prepare_staged_integration(f$theta, f$prepared, joint, free)
    expect_false(identical(joint$integration_mesh$outer[[2L]]$route, 'independent'))
    coupled = evaluate_model(f$theta, f$prepared, joint, free)
    expect_gt(max(abs(coupled$score[2:3, 'cor_re'])), 1e-5)
  }
})

test_that('native GH with staged PLN inner rules retains never-response panels and independent coordinate moments', {
  f = staged_reduction_fixture('pln')
  fixed = c(cor_re = 0, cor_error = 0); f$theta[names(fixed)] = fixed
  control = staged_test_budget(panel_control(quadrature = 'gh', nodes = c(2L, 3L),
    observation_quadrature = 'gl', observation_nodes = 12L,
    observation_options = list(reference_nodes = 15L)))
  control$integration_mesh = prepare_staged_integration(f$theta, f$prepared, control, fixed)
  actual = evaluate_model(f$theta, f$prepared, control, fixed)
  u = normal_rule(2L, 1L); v = normal_rule(3L, 1L)
  for (i in seq_len(f$prepared$n_panels)) {
    rows = which(f$prepared$group == i); selected = rows[f$prepared$d[rows] == 1]
    s = staged_selection_reference(f$theta, f$prepared, rows, staged_reduction_weights(normal = u))
    y = if (length(selected)) staged_count_reference(f$theta, f$prepared, selected,
      staged_reduction_weights(normal = v), control$integration_mesh$inner) else
      list(value = 0, score = f$theta*0)
    reference = s$score + y$score; reference[names(fixed)] = 0
    expect_equal(actual$loglik[i], s$value + y$value, tolerance = 2e-12)
    expect_equal(drop(actual$score[i, ]), reference, tolerance = 2e-11)
  }
  expect_true(all(is.finite(actual$loglik)))
  expect_equal(actual$diagnostics$panels[[1L]]$nested_error, 0)
  irrelevant = grep('^outcome\\.|^sd_re_outcome$|^sd_error$|^cor_', names(f$theta), value = TRUE)
  expect_equal(actual$score[1L, irrelevant], rep(0, length(irrelevant)), ignore_attr = TRUE)
  expected = c(sum(exp(u$log_weights)*abs(u$nodes[, 1L])),
    sum(exp(v$log_weights)*abs(v$nodes[, 1L])))
  for (i in 2:3) {
    expect_equal(actual$diagnostics$finite_prior_mass[i], 1, tolerance = 2e-14)
    expect_equal(unname(actual$diagnostics$finite_prior_absolute_moment[i, ]), expected, tolerance = 2e-14)
  }
  # Adaptive outer GH also takes its native standalone selection reduction.
  adaptive = control; adaptive$quadrature = 'agh'; adaptive$nodes = 5L
  adaptive$integration_mesh = prepare_staged_integration(f$theta, f$prepared, adaptive, fixed)
  value = evaluate_model(f$theta, f$prepared, adaptive, fixed)
  expect_true(all(is.finite(value$loglik)))
  expect_equal(value$diagnostics$panels[[1L]]$nested_error, 0)
})

test_that('independent count components preserve explicit zero-scale analytical limits', {
  for (family in c('poisson', 'pln')) for (scale in c('sd_re_selection', 'sd_re_outcome')) {
    f = staged_reduction_fixture(family)
    fixed_names = intersect(c('cor_re', 'cor_error', scale), names(f$theta))
    fixed = setNames(rep(0, length(fixed_names)), fixed_names)
    f$theta[names(fixed)] = fixed
    control = staged_reduction_control()
    control$integration_mesh = prepare_staged_integration(f$theta, f$prepared, control, fixed)
    actual = evaluate_model(f$theta, f$prepared, control, fixed)
    dimensions = c(selection = as.integer(scale != 'sd_re_selection'), outcome = as.integer(scale != 'sd_re_outcome'))
    for (i in 2:3) {
      mesh = control$integration_mesh$outer[[i]]
      expect_identical(mesh$route, 'independent')
      expect_identical(ncol(mesh$selection$lower), dimensions[['selection']])
      expect_identical(ncol(mesh$outcome$lower), dimensions[['outcome']])
      rows = which(f$prepared$group == i); selected = rows[f$prepared$d[rows] == 1]
      s = staged_selection_reference(f$theta, f$prepared, rows, staged_reduction_weights(mesh$selection))
      y = staged_count_reference(f$theta, f$prepared, selected, staged_reduction_weights(mesh$outcome),
        control$integration_mesh$inner)
      reference = s$score + y$score; reference[names(fixed)] = 0
      expect_equal(actual$loglik[i], s$value + y$value, tolerance = 2e-12)
      expect_equal(drop(actual$score[i, ]), reference, tolerance = 2e-11)
      expect_equal(actual$diagnostics$panels[[i]]$factor_dimensions, dimensions)
      expect_equal(actual$diagnostics$panels[[i]]$nodes,
        staged_node_count(mesh$selection) + staged_node_count(mesh$outcome))
    }
    expect_equal(actual$score[, scale], rep(0, f$prepared$n_panels), ignore_attr = TRUE)
  }
})
