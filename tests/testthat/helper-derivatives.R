seed_jets = function(values, q = 2L, order = 3L) {
  lapply(seq_along(values), function(j) {
    out = jet(values[j], length(values), q, order)
    out$gradient[1, j] = 1
    out
  })
}
