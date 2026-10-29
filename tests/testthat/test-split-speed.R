# Independent materialized reduction of a finite tensor. This reference remains
# in tests; production retains only the streaming generic C++ reduction.
split_bilinear_reference = function(a, b, x, y, cross, moments = TRUE) {
  terms = outer(a, b, '+') + cross*tcrossprod(x, y)
  high = max(terms); mass = exp(terms - high); total = sum(mass)
  result = list(value = high + log(total), mu = numeric(), mv = numeric(),
    weighted_v = numeric(), weighted_u = numeric())
  if (moments) {
    mass = mass/total
    result$mu = rowSums(mass); result$mv = colSums(mass)
    result$weighted_v = drop(mass %*% y)
    result$weighted_u = drop(crossprod(x, mass))
  }
  result
}

test_that('streaming tensor mass and moments equal a materialized finite sum', {
  cases = list(list(a = c(-3, -.7, -Inf), b = c(-.2, -2),
    x = c(-2, .1, 3), y = c(-1, 2), cross = -.35),
    list(a = c(-1001, -1003), b = c(-2, -4, -.3, -8),
      x = c(-4, .25), y = c(2, -1, 0, .8), cross = .7),
    list(a = c(-.3, -2, -5), b = c(-.4, -Inf),
      x = c(-3, 0, 4), y = c(-2, 1), cross = 0))
  for (z in cases) for (moments in c(FALSE, TRUE)) {
    actual = do.call(bilinear_reduce_cpp, c(z, list(moments = moments)))
    reference = do.call(split_bilinear_reference, c(z, list(moments = moments)))
    expect_equal(actual, reference, tolerance = 2e-13)
    expect_identical(actual$value, do.call(bilinear_reduce_cpp, c(z, list(moments = !moments)))$value)
  }
  expect_error(bilinear_reduce_cpp(0, 0, numeric(), 0, 0), 'dimensions')
  expect_error(bilinear_reduce_cpp(-Inf, -Inf, 0, 0, 0), 'no finite mass')
  expect_error(bilinear_reduce_cpp(NaN, 0, 0, 0, 0), 'nonfinite')
  expect_error(bilinear_reduce_cpp(0, Inf, 0, 0, 0), 'nonfinite')
  expect_error(bilinear_reduce_cpp(0, 0, 0, Inf, 0), 'nonfinite')
})

test_that('a shared partition preserves every finite node weight and cut derivative', {
  index = c(.173, -.289, .173)
  partition = split_partition(index, 1.37, 12, TRUE)
  for (order in c(3L, 12L, 15L)) {
    direct = split_axis(index, 1.37, order, 12)
    shared = split_axis(index, 1.37, order, 12, partition = partition)
    expect_identical(shared, direct)
    expect_equal(sum(exp(shared$log_weights)), 24, tolerance = 1e-13)
  }
  expect_error(split_axis(c(0, 1e-16), 1, 12L, 12), 'too close')
})

split_speed_fixture = function() {
  data = data.frame(id = rep(1:4, each = 3), d = rep(c(1, 0, 1), 4),
    y = c(0, NA, 1, 1, NA, 0, 0, NA, 1, 1, NA, 0),
    x = c(-.8, .2, .9, -.8, .2, .9, -2.4, .6, 2.7, -2.4, .6, 2.7),
    os = .043, oy = -.021)
  prepared = prepare_data(y ~ x + offset(oy), d ~ x + offset(os), data, 'id', 'probit', 'fail')
  theta = c('selection.(Intercept)' = .173, selection.x = 0,
    'outcome.(Intercept)' = -.237, outcome.x = 0,
    sd_re_selection = 1.37, sd_re_outcome = 2.13, cor_re = -.31, cor_error = 0)
  prepared$split_cache = prepare_split_panels(prepared)
  list(prepared = prepared, theta = theta)
}

test_that('within-call equation sharing retains distinct design scores and native gates', {
  z = split_speed_fixture(); cache = z$prepared$split_cache
  expect_identical(cache$keys[[1]]$selection, cache$keys[[2]]$selection)
  expect_false(identical(cache$keys[[1]]$selection_axis, cache$keys[[3]]$selection_axis))
  control = panel_control(quadrature = 'gl', nodes = 12, quadrature_options = list(reference_nodes = 15))
  actual = evaluate_split_binary(z$theta, z$prepared, control, c(cor_error = 0),
    single_outcome = FALSE)
  raw = lapply(cache$panels, function(panel) split_physical_panel(z$theta, panel, 12L,
    control$quadrature_options, fixed = c(cor_error = 0)))
  expect_equal(actual$loglik, unname(vapply(raw, `[[`, numeric(1), 'loglik')[cache$map]), tolerance = 1e-13)
  raw_score = do.call(rbind, lapply(raw, `[[`, 'score'))[cache$map, ]
  rownames(raw_score) = NULL
  expect_equal(actual$score, raw_score, tolerance = 1e-12)
  expect_gt(abs(actual$score[1, 'selection.x'] - actual$score[3, 'selection.x']), 1e-4)
  # Replace only the generic reduction, so the complete native score, tail,
  # branch and aggregate gates are compared on the identical finite rule.
  trial = new.env(parent = asNamespace('PanelSelect'))
  trial$bilinear_reduce_cpp = split_bilinear_reference
  trial$split_physical_panel = get('split_physical_panel', envir = asNamespace('PanelSelect'))
  environment(trial$split_physical_panel) = trial
  trial$evaluate_split_binary = get('evaluate_split_binary', envir = asNamespace('PanelSelect'))
  environment(trial$evaluate_split_binary) = trial
  reference = trial$evaluate_split_binary(z$theta, z$prepared, control, c(cor_error = 0),
    single_outcome = FALSE)
  expect_equal(actual$loglik, reference$loglik, tolerance = 1e-12)
  expect_equal(actual$score, reference$score, tolerance = 2e-11)
  expect_equal(actual$reference_score, reference$reference_score, tolerance = 2e-11)
  expect_equal(actual$diagnostics$aggregate_absolute_loglik_change,
    reference$diagnostics$aggregate_absolute_loglik_change, tolerance = 1e-12)
  expect_equal(actual$diagnostics$aggregate_absolute_score_change,
    reference$diagnostics$aggregate_absolute_score_change, tolerance = 2e-11)
})

test_that('zero-rho scores reuse first-pass posterior means without changing the finite rule', {
  z = split_speed_fixture(); z$theta['cor_re'] = 0
  panel = z$prepared$split_cache$panels[[1]]
  options = list(bound = 12, max_nodes = 1000000L)
  mean_node = function(index, signs, scale) {
    axis = split_axis(index, scale, 12L, 12, FALSE)
    lp = split_probit_values(index, signs, scale, axis$nodes, FALSE)$log_product
    terms = axis$log_weights + dnorm(axis$nodes, log = TRUE) + lp
    sum(axis$nodes*exp(terms - split_logsum(terms)))
  }
  a = drop(panel$w %*% z$theta[1:2]) + panel$offset_selection
  b = drop(panel$x %*% z$theta[3:4]) + panel$offset_outcome
  value = split_physical_panel(z$theta, panel, 12L, options, fixed = c(cor_error = 0))
  expected = mean_node(a, panel$selection_sign, z$theta['sd_re_selection'])*
    mean_node(b, panel$outcome_sign, z$theta['sd_re_outcome'])
  expect_identical(unname(value$score['cor_re']), expected)
  expect_gt(abs(expected), 1e-5)
})

test_that('private reuse binds equations, layouts, fixed masks, scores and changed prepared inputs', {
  z = split_speed_fixture(); prepared = z$prepared
  prepared$split_workspace = new_split_workspace(prepared)
  workspace = prepared$split_workspace
  control = panel_control(quadrature = 'gl', nodes = 12, quadrature_options = list(reference_nodes = 15))
  run = function(theta = z$theta, data = prepared, fixed = c(cor_error = 0), score = TRUE)
    evaluate_split_binary(theta, data, control, fixed, score, single_outcome = FALSE)
  cold = run(); context = workspace$context
  warm = run(); expect_identical(warm$score, cold$score)
  expect_identical(workspace$context$axes, context$axes)
  selection_entries = grep('^selection:', ls(context$axes), value = TRUE)
  kept = get(selection_entries[1], context$axes)
  changed = z$theta; changed['outcome.(Intercept)'] = changed['outcome.(Intercept)'] + .113
  different = run(changed)
  expect_identical(get(selection_entries[1], context$axes), kept)
  fresh = prepared; fresh$split_workspace = NULL
  expect_equal(different$score, run(changed, fresh)$score, tolerance = 1e-12)
  for (change in c('offset', 'design')) {
    altered = prepared
    if (change == 'offset') altered$offset_selection[1] = altered$offset_selection[1] + .087 else
      altered$w[1, 'x'] = altered$w[1, 'x'] + .123
    shared = run(data = altered)
    fresh = altered; fresh$split_workspace = fresh$split_cache = NULL
    expect_equal(shared$loglik, run(data = fresh)$loglik, tolerance = 1e-12)
    expect_equal(shared$score, run(data = fresh)$score, tolerance = 1e-12)
  }
  fixed = c(cor_error = 0, selection.x = 0)
  expect_equal(run(fixed = fixed)$score[, 'selection.x'], rep(0, 4))
  expect_null(run(score = FALSE)$score)
  expect_identical(run(score = FALSE)$loglik, run()$loglik)
  reversed = z$theta[rev(seq_along(z$theta))]
  expect_equal(run(reversed)$score[, names(z$theta)], run()$score, tolerance = 1e-12)
})

test_that('split prepared signatures retain counts and factorials and bound cache storage', {
  data = data.frame(id = rep(1:3, each = 2), y = c(0, 1, 0, 2, 0, 1))
  for (family in c('poisson', 'pln')) {
    cache = prepare_split_panels(prepare_data(y ~ 1, NULL, data, 'id', family))
    expect_identical(cache$map, c(1L, 2L, 1L))
    expect_identical(cache$model_id, paste0(family, '_re'))
    expect_identical(cache$panels[[2]]$outcome, c(0, 2))
    expect_equal(cache$panels[[2]]$log_factorial, lgamma(c(0, 2) + 1))
  }
  memo = new.env(parent = emptyenv())
  for (j in seq_len(70L)) split_reuse(memo, as.character(j), function() j)
  expect_lte(length(memo), 64L)
  expect_lte(attr(memo, 'bytes'), 16*1024^2)
  expect_identical(split_reuse(memo, '70', function() stop('must reuse')), 70L)
})

test_that('private Gaussian reuse binds changed outcomes, offsets and derivative designs', {
  data = data.frame(id = rep(1:2, each = 3), d = rep(c(1, 0, 1), 2),
    y = c(-.7, NA, .3, .4, NA, -.2), x = rep(c(-.83, .27, .91), 2), oy = .051)
  prepared = prepare_data(y ~ x + offset(oy), d ~ x, data, 'id', 'gaussian', 'fail')
  prepared$split_cache = prepare_split_panels(prepared)
  prepared$split_workspace = new_split_workspace(prepared)
  theta = c('selection.(Intercept)' = .173, selection.x = -.283,
    'outcome.(Intercept)' = -.237, outcome.x = .193,
    sd_re_selection = 1.37, sd_re_outcome = 2.13, sd_error = .9, cor_re = -.31, cor_error = 0)
  control = panel_control(quadrature = 'gl', nodes = 12, quadrature_options = list(reference_nodes = 15))
  run = function(at = theta, p = prepared)
    evaluate_split_binary(at, p, control, c(cor_error = 0))
  first = run(); marginal = prepared$split_workspace$gaussian_value
  different = theta; different['selection.(Intercept)'] = different['selection.(Intercept)'] + .117
  run(different)
  expect_identical(prepared$split_workspace$gaussian_value, marginal)
  for (change in c('outcome', 'offset', 'design')) {
    altered = prepared
    if (change == 'outcome') altered$y[1] = altered$y[1] + .121 else
      if (change == 'offset') altered$offset_outcome[1] = altered$offset_outcome[1] + .073 else
        altered$x[1, 'x'] = altered$x[1, 'x'] + .087
    shared = run(p = altered)
    fresh = altered
    fresh$model_data = fresh$integration_cache = fresh$split_cache = fresh$split_workspace = NULL
    expected = run(p = fresh)
    expect_equal(shared$loglik, expected$loglik, tolerance = 1e-12)
    expect_equal(shared$score, expected$score, tolerance = 1e-12)
    expect_gt(max(abs(shared$loglik - first$loglik)), 1e-5)
  }
})
