# The shared assembly takes a panel's row sum either from group_sum() or from
# the final (weight) column of group_weighted_design(). Both sum each panel's
# rows in row order with Armadillo's two-accumulator reduction. Finite sums,
# signed zeros and infinities agree bitwise. Invalid sums can propagate either
# NA or NaN when both occur, so compare their missingness rather than their encoding.
test_that('panel sums from group_sum and the design weight column agree', {
  local_test_seed(29)
  special = c(0, -0, 4.9e-324, -4.9e-324, 1e-310, 1, -1, 1/3, 1e300, -1e300, Inf, -Inf, NaN, NA)
  for (count in c(0L, 1L, 2L, 3L, 4L, 5L, 7L, 100L)) for (blocks in c(1L, 3L)) for (kind in 1:3) {
    # Four panels, some of them empty.
    lengths = tabulate(sample.int(4L, count, replace = TRUE), nbins = 4L)
    ptr = c(0L, cumsum(lengths))
    size = count*blocks
    values = switch(kind, stats::rnorm(size), sample(special, size, replace = TRUE),
      stats::rnorm(size)*10^sample(-300:300, size, replace = TRUE))
    x = matrix(values, count, blocks)
    design = matrix(stats::rnorm(count*2L), count, 2L)
    product = group_weighted_design(design, x, ptr)
    actual = product[, ncol(product)]; expected = as.vector(group_sum(x, ptr))
    label = paste('rows', count, 'blocks', blocks, 'values', kind)
    expect_identical(is.na(actual), is.na(expected), label = label)
    expect_true(identical(actual[!is.na(actual)], expected[!is.na(expected)], num.eq = FALSE),
      info = label)
  }
})
