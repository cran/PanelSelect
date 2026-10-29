# One rule for length-two nodes: element 1 orders the selection effect and
# element 2 the outcome effect on every route, including when one effect remains.

nodes_fixture = function(seed = 2201) {
  local_test_seed(seed)
  panels = 12L; waves = 3L; n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels)
  data$d = as.integer(.2 + .6*data$w + .7*u[data$id] + stats::rnorm(n) > 0)
  data$y = as.integer(.1 + .5*data$x + .6*u[data$id] + stats::rnorm(n) > 0)
  data$count = stats::rpois(n, exp(.2 + .3*data$x + .5*u[data$id]))
  data$y[data$d == 0] = NA
  data
}

nodes_theta = c('selection.(Intercept)' = .2, selection.w = .6, 'outcome.(Intercept)' = .1,
  outcome.x = .5, sd_re_selection = .8, sd_re_outcome = .7, cor_re = .3, cor_error = .3)

evaluate_nodes = function(data, outcome, selection, family, fixed, nodes, quadrature = 'gh') {
  prepared = prepare_data(outcome, selection, data, 'id', family)
  info = parameter_layout(prepared)
  fixed = resolve_fixed(fixed, info, prepared)
  theta = stats::setNames(rep(.25, nrow(info)), info$name)
  shared = intersect(names(nodes_theta), info$name)
  theta[shared] = nodes_theta[shared]
  theta[names(fixed)] = fixed
  evaluate_model(theta, prepared, panel_control(quadrature = quadrature, nodes = nodes), fixed)
}

test_that('a remaining outcome effect uses the second order on joint and history routes', {
  data = nodes_fixture()
  fixed = c(sd_re_selection = 0)
  for (quadrature in c('gh', 'agh')) {
    run = function(nodes) evaluate_nodes(data, y ~ x, d ~ w, 'probit', fixed, nodes, quadrature)
    pair = run(c(3L, 7L)); seven = run(7L); three = run(3L)
    expect_identical(pair$loglik, seven$loglik)
    expect_identical(pair$score, seven$score)
    expect_false(isTRUE(all.equal(pair$loglik, three$loglik)))
    expect_true(all(quadrature_outer_orders(pair$diagnostics) == 7L))
  }
  theta = nodes_theta; theta[c('sd_re_selection', 'cor_re')] = 0
  fit = function(nodes) panel_select(y ~ x, d ~ w, data, 'id', 'probit', fixed = theta,
    control = panel_control(quadrature = 'gh', nodes = nodes))
  targets = data[data$id %in% c(2, 5), ]
  history = function(object) predict(object, targets, conditioning = 'history', history = data)
  expect_identical(history(fit(c(3L, 7L))), history(fit(7L)))
})

test_that('standalone and selection-only effects take their own orders', {
  data = nodes_fixture()
  standalone = function(nodes) evaluate_nodes(data, count ~ x, NULL, 'poisson', NULL, nodes)
  expect_identical(standalone(c(3L, 7L))$loglik, standalone(7L)$loglik)
  selection = function(nodes) evaluate_nodes(data, y ~ x, d ~ w, 'probit', c(sd_re_outcome = 0), nodes)
  expect_identical(selection(c(3L, 7L))$loglik, selection(3L)$loglik)
  both = function(nodes) evaluate_nodes(data, y ~ x, d ~ w, 'probit', NULL, nodes)
  expect_false(isTRUE(all.equal(both(c(3L, 7L))$loglik, both(c(7L, 3L))$loglik)))
})

test_that('named nodes are stored in canonical order and partial names fail', {
  expect_identical(panel_control(nodes = c(outcome = 7, selection = 3))$nodes, c(3L, 7L))
  expect_identical(panel_control(nodes = c(selection = 3, outcome = 7)), panel_control(nodes = c(3, 7)))
  expect_error(panel_control(nodes = c(sel = 3, outcome = 7)), 'Named nodes must be')
  expect_error(panel_control(nodes = c(selection = 3)), 'Named nodes must be')
  data = nodes_fixture()
  named = evaluate_nodes(data, y ~ x, d ~ w, 'probit', c(sd_re_selection = 0),
    c(outcome = 7L, selection = 3L))
  expect_identical(named$loglik, evaluate_nodes(data, y ~ x, d ~ w, 'probit',
    c(sd_re_selection = 0), c(3L, 7L))$loglik)
})
