# Independently enumerate positive two-/three-row products. These tests verify
# algebraic perturbation/tail envelopes, not numerical integration accuracy.
staged_product_fixture = function(rows, outer_order = 2L) {
  data = data.frame(id = 1L, y = c(0, 1, 4)[seq_len(rows)],
    x = c(-.4, .7, 1.1)[seq_len(rows)])
  prepared = prepare_data(y ~ x, NULL, data, 'id', 'pln', 'fail')
  theta = c('outcome.(Intercept)' = -.237, outcome.x = .193,
    sd_re_outcome = .61, sd_error = .53)
  control = panel_control(quadrature = 'gh', nodes = 2L,
    observation_quadrature = 'gh', observation_nodes = 3L)
  context = make_model_context(theta, prepared, control, numeric())
  bound = c(4, 5, 6)[seq_len(rows)]
  context$inner$quadrature = 'gl'
  context$inner$meshes = lapply(bound, function(b) list(bound = b))
  # Explicit standard-normal GH2/GH3 rules. Their first absolute moments are
  # 1 and sqrt(3)/3, respectively, and need not equal sqrt(2/pi).
  outer = if (outer_order == 2L) list(z = c(-1, 1), weight = c(.5, .5)) else
    list(z = c(-sqrt(3), 0, sqrt(3)), weight = c(1/6, 2/3, 1/6))
  outer$weight = 1.17*outer$weight
  inner_mass = c(1.25, 1.07, 1.60)[seq_len(rows)]
  actual_first = c(1.25, 1.07*sqrt(3)/3, 1.60*sqrt(3)/3)[seq_len(rows)]
  first = pmax(sqrt(2/pi), actual_first)
  y = prepared$y; maximum = dpois(y, y)
  derivative = y*dpois(y, y)+(y+1)*dpois(y+1, y+1)
  C = maximum*pmax(1, inner_mass)
  n = length(outer$z); p = length(theta)
  truth = finite = matrix(0, n, rows)
  dtruth = dfinite = array(0, c(n, rows, p))
  for (j in seq_len(n)) for (t in seq_len(rows)) {
    truth[j, t] = maximum[t]*(.25+.12*((j+t) %% 4))
    finite[j, t] = C[t]*(.21+.11*((2*j+t) %% 5))
    D = c(derivative[t]*max(1, inner_mass[t]),
      abs(data$x[t])*derivative[t]*max(1, inner_mass[t]),
      abs(outer$z[j])*derivative[t]*max(1, inner_mass[t]),
      derivative[t]*first[t])
    dtruth[j, t, ] = D*sin(j*.37+t*.29+seq_len(p)*.17)
    dfinite[j, t, ] = D*cos(j*.23+t*.41+seq_len(p)*.13)
  }
  product = function(value, gradient) {
    mass = apply(value, 1L, prod)
    first = matrix(0, n, p)
    for (j in seq_len(n)) for (t in seq_len(rows))
      first[j, ] = first[j, ]+gradient[j, t, ]*prod(value[j, -t])
    list(value = sum(outer$weight*mass),
      first = colSums(outer$weight*first))
  }
  old = product(truth, dtruth); new = product(finite, dfinite)
  error = apply(abs(finite-truth), 2L, max)
  score_error = matrix(0, rows, p)
  for (t in seq_len(rows)) for (k in seq_len(p))
    score_error[t, k] = max(abs(dfinite[, t, k]-dtruth[, t, k]))
  recorder = new.env(parent = emptyenv())
  recorder$outer_mass = sum(outer$weight)
  recorder$outer_moment = matrix(c(sum(outer$weight*abs(outer$z)), 0), 1L, 2L)
  recorder$inner_mass = inner_mass; recorder$inner_first = first
  recorder$inner_error = log(error); recorder$inner_score_error = log(score_error)
  recorder$inner_probes = recorder$inner_score_probes = vector('list', rows)
  list(theta = theta, prepared = prepared, context = context, recorder = recorder,
    truth = truth, finite = finite, old = old, new = new, product = product,
    dtruth = dtruth, dfinite = dfinite, C = C, error = error,
    score_error = score_error, maximum = maximum, bound = bound, outer = outer)
}

test_that('count-aware nested envelopes bound independently enumerated positive products', {
  for (rows in 2:3) for (outer_order in 2:3) {
    f = staged_product_fixture(rows, outer_order)
    score = f$new$first/f$new$value
    actual = staged_nested_error(f$theta, f$prepared, f$context,
      log(f$new$value), score, f$recorder)
    companion = vapply(seq_len(rows), function(t) prod(f$C[-t]), numeric(1L))
    E = sum(f$error*companion)
    expected_probability = sum(f$outer$weight)*E/f$new$value
    expect_equal(actual$probability, expected_probability, tolerance = 1e-12)
    expect_gte(actual$probability*f$new$value,
      abs(f$new$value-f$old$value)-1e-14)
    # The direct numerator envelope is separate from the normalization term.
    numerator = (actual$score-abs(score)*actual$probability)*f$new$value
    expect_true(all(numerator+1e-12 >= abs(f$new$first-f$old$first)))
    # Companion derivative terms must remain present even with identical local
    # derivative numerators. Only probabilities change in this second fixture.
    f$recorder$inner_score_error[,] = -Inf
    same_first = f$product(f$finite, f$dtruth)
    second = staged_nested_error(f$theta, f$prepared, f$context,
      log(same_first$value), same_first$first/same_first$value, f$recorder)
    numerator2 = (second$score-abs(same_first$first/same_first$value)*
      second$probability)*same_first$value
    expect_true(all(numerator2+1e-12 >= abs(same_first$first-f$old$first)))
    expect_true(any(f$recorder$inner_mass > 1))
    expect_gt(f$recorder$inner_first[1L], sqrt(2/pi))
  }
})

test_that('analytic observation tails receive the same count companion bounds', {
  for (rows in 2:3) {
    f = staged_product_fixture(rows, 3L)
    actual = staged_nested_error(f$theta, f$prepared, f$context,
      log(f$new$value), f$new$first/f$new$value, f$recorder)
    tail_mass = 2*pnorm(-f$bound)
    companion = vapply(seq_len(rows), function(t) prod(f$C[-t]), numeric(1L))
    expected = sum(f$outer$weight)*sum(f$maximum*tail_mass*companion)/f$new$value
    expect_equal(actual$tail, expected, tolerance = 1e-12)
    # Enumerate an admissible positive omitted tail of each exact row. The
    # enlarged true rows still lie under their maxima and hence under C_t.
    enlarged = sweep(f$truth, 2L, .25*f$maximum*tail_mass, '+')
    omitted = f$product(enlarged, f$dtruth)$value-f$old$value
    expect_gte(actual$tail*f$new$value, omitted-1e-14)
    expect_true(all(enlarged <= rep(f$C, each = nrow(enlarged))))
  }
})
