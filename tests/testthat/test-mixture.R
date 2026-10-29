test_that('batched log mixtures agree with sequential mixtures through order three', {
  for (order in 0:3) {
    variables = seed_jets(c(0.2, -0.3, 0.4), order = order)
    n = 7L; nodes = c(-2, -0.5, 0.4, 2)
    row = rep(seq_len(n), length(nodes))
    node = rep(nodes, each = n)
    eta = variables[[1]] + variables[[2]]*row/10
    terms = -.5*((node-eta)/exp(variables[[3]]))^2 - variables[[3]]
    sequential = jet(rep(-Inf, n), terms$d, terms$q, terms$order)
    for (k in seq_along(nodes)) sequential = jet_logadd(sequential,
      jet_subset(terms, seq_len(n)+(k-1L)*n))
    batch = jet_log_mixture(terms, n)
    expect_equal(batch, sequential, tolerance = 1e-11)
    blocked = jet(rep(0, n), terms$d, terms$q, terms$order)
    for (rows in list(1:2, 3:5, 6:7)) {
      index = unlist(lapply(seq_along(nodes), function(k) rows+(k-1L)*n))
      blocked = jet_replace(blocked, rows, jet_log_mixture(jet_subset(terms, index), length(rows)))
    }
    expect_equal(blocked, batch, tolerance = 1e-12)
    shifted = terms + rep(seq(-1000, 1000, length.out = n), length(nodes))
    stabilized = jet_log_mixture(shifted, n)
    expect_equal(stabilized$value, batch$value+seq(-1000, 1000, length.out = n), tolerance = 1e-12)
    expect_equal(stabilized$gradient, batch$gradient, tolerance = 1e-10)
    if (order >= 3L) expect_equal(stabilized$third, batch$third, tolerance = 1e-10)
  }
})

test_that('zero masses do not introduce infinite log derivatives', {
  term = jet(c(-Inf, -Inf, 0.3, -Inf), 1L, 1L, 3L)
  term$gradient[, 1] = c(Inf, Inf, 2, Inf)
  term$hessian[, 1, 1] = c(Inf, Inf, 3, Inf)
  term$third[, 1, 1, 1] = c(Inf, Inf, 4, Inf)
  out = jet_log_mixture(term, 2L)
  expect_equal(out$value, c(0.3, -Inf))
  expect_equal(out$gradient[, 1], c(2, 0))
  expect_equal(out$hessian[, 1, 1], c(3, 0))
  expect_equal(out$third[, 1, 1, 1], c(4, 0))
})

test_that('invalid mixture components cannot be silently omitted', {
  term = jet(c(NaN, Inf, 1000, 0), 1L, 1L, 3L)
  out = jet_log_mixture(term, 2L)
  expect_true(all(!is.finite(out$value)))
  expect_true(all(!is.finite(out$gradient)))
})

test_that('mixture batches partition observations and bound larger derivative work', {
  context = jet(0, 200L, 2L, 3L)
  groups = mixture_batches(1000L, 31L, context)
  expect_gt(length(groups), 1L)
  expect_identical(unname(unlist(groups)), seq_len(1000L))
  expect_lt(max(lengths(groups)), 1000L)
  expect_length(mixture_batches(0L, 31L, context), 0L)
})
