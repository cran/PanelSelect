physical_pln_fixture = function() {
  data = data.frame(id = rep('one', 3L), d = c(1,0,1), y = c(1,NA,3),
    x = c(-.57,.41,1.13), w = c(.37,-.59,.73),
    os = c(-.057,.023,.089), oy = c(-.043,.033,.067))
  prepared = prepare_data(y ~ x + offset(oy), d ~ w + offset(os),
    data, 'id', 'pln', 'fail')
  theta = c('selection.(Intercept)' = .173, selection.w = -.283,
    'outcome.(Intercept)' = -.237, outcome.x = .193,
    sd_re_selection = 1.25, sd_re_outcome = 1.7, sd_error = .65, cor_re = -.3, cor_error = 0)
  list(theta = theta[parameter_layout(prepared)$name], prepared = prepared)
}

physical_pln_control = function() staged_test_budget(panel_control(quadrature = 'gl',
  nodes = 15L, observation_quadrature = 'gl', observation_nodes = 15L,
  quadrature_options = list(reference_nodes = 20L),
  observation_options = list(reference_nodes = 20L)))

# Independent finite product and derivative calculation. No production model
# kernels, compiled bilinear reducers or physical score assembler are called.
physical_pln_finite_reference = function(theta, prepared, control, fixed, reference = FALSE) {
  mesh = control$integration_mesh$outer[[1L]]
  u = staged_cell_rule(mesh$u, reference = reference)
  v = staged_cell_rule(mesh$v, reference = reference)
  U = u$nodes[, 1L]; V = v$nodes[, 1L]; rho = theta[['cor_re']]; q = 1-rho^2
  s_names = paste0('selection.', colnames(prepared$w)); y_names = paste0('outcome.', colnames(prepared$x))
  a = drop(prepared$w %*% theta[s_names])+prepared$offset_selection
  b = drop(prepared$x %*% theta[y_names])+prepared$offset_outcome
  signs = 2*prepared$d-1
  index = sweep(outer(U, rep(theta[['sd_re_selection']], length(a))), 2L, a, '+')
  signed = sweep(index, 2L, signs, '*')
  logp = matrix(pnorm(signed, log.p = TRUE), length(U), length(a))
  part = sweep(exp(dnorm(signed, log = TRUE)-logp), 2L, signs, '*')
  selection = rowSums(logp)
  sg = matrix(0, length(U), length(theta), dimnames = list(NULL,names(theta)))
  sg[, s_names] = part %*% prepared$w; sg[, 'sd_re_selection'] = U*rowSums(part)
  outcome = numeric(length(V)); yg = matrix(0, length(V), length(theta),
    dimnames = list(NULL,names(theta)))
  for (row in which(prepared$d == 1)) {
    rule = staged_inner_rule(control$integration_mesh$inner[[row]], reference)
    E = rule$nodes[, 1L]
    eta = outer(theta[['sd_re_outcome']]*V+b[row], theta[['sd_error']]*E, '+')
    mean = exp(eta)
    log_terms = matrix(dpois(prepared$y[row], mean, log = TRUE), nrow(mean), ncol(mean))+
      rep(rule$log_weights, each = length(V))
    high = apply(log_terms, 1L, max); weights = exp(log_terms-high)
    total = rowSums(weights); weights = weights/total
    outcome = outcome+high+log(total)
    db = rowSums(weights*(prepared$y[row]-mean))
    de = rowSums(sweep(weights*(prepared$y[row]-mean), 2L, E, '*'))
    yg[, y_names] = yg[, y_names, drop = FALSE]+outer(db, prepared$x[row, ])
    yg[, 'sd_re_outcome'] = yg[, 'sd_re_outcome']+V*db
    yg[, 'sd_error'] = yg[, 'sd_error']+de
  }
  log_terms = outer(u$log_weights+selection-U^2/(2*q),
    v$log_weights+outcome-V^2/(2*q), '+')+outer(U,(rho/q)*V)-log(2*pi)-.5*log(q)
  high = max(log_terms); mass = exp(log_terms-high); total = sum(mass); mass = mass/total
  score = drop(crossprod(rowSums(mass), sg)+crossprod(colSums(mass), yg))
  names(score) = names(theta)
  score['cor_re'] = rho/q+sum(mass*((1+rho^2)*outer(U,V)-
    rho*outer(U^2,V^2,'+')))/q^2
  score[names(fixed)] = 0
  list(loglik = high+log(total), score = score)
}

test_that('physical PLN reuses immutable meshes with exact canonical and reference finite scores', {
  fixture = physical_pln_fixture(); fixed = c(cor_error = 0)
  control = physical_pln_control()
  control$integration_mesh = prepare_staged_integration(fixture$theta, fixture$prepared, control, fixed)
  expect_identical(control$integration_mesh$outer[[1L]]$route, 'physical_pln')
  saved = unserialize(serialize(control$integration_mesh, NULL))
  actual = evaluate_model(fixture$theta, fixture$prepared, control, fixed)
  for (reference in c(FALSE,TRUE)) {
    expected = physical_pln_finite_reference(fixture$theta, fixture$prepared, control, fixed, reference)
    if (reference) expect_equal(drop(actual$reference_score), expected$score, tolerance = 2e-11) else {
      expect_equal(actual$loglik, expected$loglik, tolerance = 2e-12)
      expect_equal(drop(actual$score), expected$score, tolerance = 2e-11)
    }
  }
  frozen = control; frozen$.staged_preparing = TRUE
  free = !names(fixture$theta) %in% names(fixed)
  numerical = core_difference(function(at) {
    theta = fixture$theta; theta[names(at)] = at
    sum(evaluate_model(theta, fixture$prepared, frozen, fixed, score = FALSE)$loglik)
  }, fixture$theta[free])
  expect_lt(max(abs(colSums(actual$score)[free]-numerical)), 2e-6)
  expect_identical(evaluate_model(fixture$theta, fixture$prepared, control, fixed,
    score = FALSE)$loglik, actual$loglik)
  likelihood_only = evaluate_model(fixture$theta, fixture$prepared, control, fixed, score = FALSE)
  likelihood_panel = likelihood_only$diagnostics$panels[[1L]]
  expect_false(likelihood_panel$score_checks_available)
  expect_false(likelihood_panel$physical_prior$score_moments_evaluated)
  expect_lt(likelihood_panel$prior_rule_node_work, actual$diagnostics$panels[[1L]]$prior_rule_node_work)
  expect_equal(likelihood_panel$tail, actual$diagnostics$panels[[1L]]$tail, tolerance = 2e-13)
  expect_equal(evaluate_model(fixture$theta, fixture$prepared, control, fixed)$score, actual$score)
  replay = unserialize(serialize(control, NULL))
  expect_equal(evaluate_model(fixture$theta, fixture$prepared, replay, fixed)$score, actual$score)
  # Compiled access may materialize deferred row names without changing the
  # mesh. Compare every value/attribute to an independent saved object; raw
  # serialization bytes also encode R's transient ALTREP representation.
  expect_identical(control$integration_mesh, saved)
  limited = control
  limited$quadrature_options$max_nodes = staged_node_count(control$integration_mesh$outer[[1L]])-1
  limited$integration_mesh$binding = staged_binding(fixture$prepared, limited, fixed)
  expect_error(evaluate_model(fixture$theta, fixture$prepared, limited, fixed),
    'node or cell budget', class = 'panelselect_integration_error')
})

test_that('physical PLN retains free-rho-zero scores and exact fixed-zero reductions', {
  fixture = physical_pln_fixture(); fixture$theta['cor_re'] = 0
  fixed = c(cor_error = 0); control = physical_pln_control()
  control$integration_mesh = prepare_staged_integration(fixture$theta, fixture$prepared, control, fixed)
  actual = evaluate_model(fixture$theta, fixture$prepared, control, fixed)
  expected = physical_pln_finite_reference(fixture$theta, fixture$prepared, control, fixed)
  expect_gt(abs(actual$score[1L,'cor_re']), 1e-5)
  expect_equal(unname(actual$score[1L,'cor_re']), expected$score[['cor_re']], tolerance = 2e-12)
  fixed = c(cor_re = 0, cor_error = 0)
  independent = physical_pln_control()
  independent$integration_mesh = prepare_staged_integration(fixture$theta, fixture$prepared, independent, fixed)
  expect_identical(independent$integration_mesh$outer[[1L]]$route, 'independent')
  expect_equal(evaluate_model(fixture$theta, fixture$prepared, independent, fixed)$loglik,
    actual$loglik, tolerance = 1e-10)
  # The occasion correlation remains estimable when free atzero.
  fixture$theta['cor_re'] = .3
  coupled = physical_pln_control(); coupled$.staged_seed = TRUE
  # This check concerns dispatch, not refinement of a costly coupled fixture.
  # Independent multiwave/extended tests cover the complete free-tau rule.
  layout = staged_panel_layout(fixture$theta, fixture$prepared, coupled, numeric())
  expect_identical(layout$route, 'joint')
  expect_false(staged_physical_pln_eligible(fixture$theta, fixture$prepared, coupled, numeric()))
  context = make_model_context(fixture$theta, fixture$prepared, coupled, numeric())
  seed = staged_seed_mesh(context, fixture$prepared, fixture$theta, coupled, numeric())
  expect_identical(ncol(seed$lower), 2L)
})

test_that('physical PLN panel mappings preserve ragged histories, offsets and row permutations', {
  data = data.frame(id = c('a','a','a','b','b','c'), d = c(1,0,1,0,0,1),
    y = c(1,NA,3,NA,NA,2), x = c(-.57,.41,1.13,.2,-.3,.4),
    w = c(.37,-.59,.73,-.1,.7,.6), os = c(-.057,.023,.089,.02,-.03,.04),
    oy = c(-.043,.033,.067,-.01,.02,.03))
  fixture = physical_pln_fixture(); theta = fixture$theta; fixed = c(cor_error = 0)
  evaluate = function(data) {
    prepared = prepare_data(y ~ x+offset(oy), d ~ w+offset(os), data, 'id', 'pln', 'fail')
    control = physical_pln_control()
    control$integration_mesh = prepare_staged_integration(theta, prepared, control, fixed)
    list(value = evaluate_model(theta, prepared, control, fixed),
      prepared = prepared, control = control)
  }
  first = evaluate(data); second = evaluate(data[c(6,4,1,5,3,2), ])
  at = match(first$prepared$id_values, second$prepared$id_values)
  expect_equal(second$value$loglik[at], first$value$loglik, tolerance = 2e-12)
  expect_equal(second$value$score[at, , drop = FALSE], first$value$score, tolerance = 2e-11)
  never = match('b', first$prepared$id_values)
  expect_identical(first$control$integration_mesh$outer[[never]]$route, 'selection_only')
  expect_equal(unname(first$value$score[never, c('sd_re_outcome','sd_error','cor_re','cor_error')]), rep(0,4))
  changed = first$prepared; changed$y[which(changed$d == 1)[1L]] = 5
  expect_error(evaluate_model(theta, changed, first$control, fixed), 'not bound')
})

test_that('physical refinement uses Cartesian cell evidence and retains positive axis rules', {
  axis = function(cuts) c(staged_cells(list(cuts)),
    list(method = 'gl', rule = staged_rule(15L, 20L, q = 1L), bound = max(abs(cuts))))
  mesh = list(route = 'physical_pln', u = axis(c(-2,0,2)), v = axis(c(-1,0,1)))
  value = list(cell_refinement_error = c(0,3,0,0), dimension_error = matrix(0,4L,2L))
  value$dimension_error[2L,2L] = 3
  refined = staged_physical_refine(mesh, value)
  expect_identical(nrow(refined$u$lower), 2L)
  expect_identical(nrow(refined$v$lower), 3L)
  expect_equal(sum(exp(staged_cell_rule(refined$u)$log_weights))*
    sum(exp(staged_cell_rule(refined$v)$log_weights)), 8, tolerance = 2e-13)
  expect_identical(staged_interval_count(refined), 6L)
  # Each cell pair holds the canonical and the reference tensor rules.
  expect_equal(staged_node_count(refined), 6*(15^2+20^2))
})
