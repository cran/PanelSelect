# Retained derivative-engine mixtures for reference verification, outside
# production likelihood dispatch. Columns are nodes; rows are observation integrals.
# The centering values are constants for differentiation: they cancel exactly
# from log(sum(exp(log_term - center))) + center.
jet_log_mixture = function(log_term, n) {
  if (log_term$order >= 2L) return(jet_log_mixture_r(log_term, n))
  score = if (log_term$order) log_term$gradient else matrix(numeric(), length(log_term$value), 0L)
  reduced = logsum_score_cpp(log_term$value, score, n)
  out = jet(reduced[, 1L], log_term$d, log_term$q, log_term$order)
  if (log_term$order) out$gradient = reduced[, -1L, drop = FALSE]
  out
}

# The R expression is also the reference for the first-score native reduction.
jet_log_mixture_r = function(log_term, n) {
  values = matrix(log_term$value, nrow = n)
  group = rep(seq_len(n), ncol(values))
  invalid = which(rowSums(is.na(values) | values == Inf) > 0)
  center = values[cbind(seq_len(n), max.col(values, ties.method = "first"))]
  center[!is.finite(center)] = 0
  mass = jet(rep(0, length(group)), log_term$d, log_term$q, log_term$order)
  finite = which(is.finite(log_term$value) & !group %in% invalid)
  if (length(finite)) mass = jet_replace(mass, finite,
    exp(jet_subset(log_term, finite) - center[group[finite]]))
  total = jet_group_sum(mass, group, n)
  result = jet(rep(-Inf, n), log_term$d, log_term$q, log_term$order)
  positive = which(total$value > 0)
  if (length(positive)) result = jet_replace(result, positive,
    log(jet_subset(total, positive)) + center[positive])
  if (length(invalid)) {
    bad = jet(rep(NaN, length(invalid)), log_term$d, log_term$q, log_term$order)
    for (field in c("gradient", "hessian", "third"))
      if (!is.null(bad[[field]])) bad[[field]][] = NaN
    result = jet_replace(result, invalid, bad)
  }
  result
}

# Limit expanded inner-node storage to roughly two million derivative entries
# per block. No observation x outer-node x parameter tensor is retained.
mixture_batches = function(n, nodes, context) {
  components = 1 + context$d*(context$order > 0) +
    context$q*context$d*(context$order > 1) +
    context$q^2*context$d*(context$order > 2)
  size = max(1L, floor(2e6/(nodes*components)))
  split(seq_len(n), ceiling(seq_len(n)/size))
}
