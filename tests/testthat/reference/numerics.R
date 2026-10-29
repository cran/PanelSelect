# Small independent R references support parity tests and paired benchmarks.
group_sum_r = function(x, ptr) {
  out = matrix(0, length(ptr)-1L, ncol(x))
  for (g in seq_len(nrow(out))) {
    rows = seq_len(ptr[g+1L]-ptr[g]) + ptr[g]
    out[g, ] = colSums(x[rows, , drop = FALSE])
  }
  out
}

group_weighted_design_r = function(x, weights, ptr, mixture_weights = NULL) {
  design = cbind(x, 1)
  groups = length(ptr)-1L
  if (!is.null(mixture_weights)) {
    collapsed = numeric(nrow(x))
    for (g in seq_len(groups)) {
      rows = seq_len(ptr[g+1L]-ptr[g])+ptr[g]
      for (b in seq_len(ncol(weights))) if (is.na(mixture_weights[g, b]) || mixture_weights[g, b] != 0)
        collapsed[rows] = collapsed[rows]+mixture_weights[g, b]*weights[rows, b]
    }
    weights = matrix(collapsed, ncol = 1L)
  }
  out = matrix(0, groups*ncol(weights), ncol(design))
  for (b in seq_len(ncol(weights))) out[seq_len(groups)+groups*(b-1L), ] =
    group_sum_r(design*weights[, b], ptr)
  out
}

mixture_weights_r = function(x) {
  if (anyNA(x) || any(x == Inf)) return(list(value = NaN, weights = rep(NaN, length(x))))
  if (!length(x) || all(x == -Inf)) return(list(value = -Inf, weights = rep(0, length(x))))
  center = max(x)
  mass = exp(x-center)
  list(value = center+log(sum(mass)), weights = mass/sum(mass))
}

logsum_weights_r = function(log_term, n) {
  nodes = if (n) length(log_term)/n else 0L
  out = list(value = rep(-Inf, n), weights = matrix(0, n, nodes))
  for (i in seq_len(n)) {
    rule = mixture_weights_r(log_term[i+n*(seq_len(nodes)-1L)])
    out$value[i] = rule$value
    out$weights[i, ] = rule$weights
  }
  out
}

logsum_score_r = function(log_term, score, n) {
  out = matrix(0, n, ncol(score)+1L)
  nodes = if (n) length(log_term)/n else 0L
  for (i in seq_len(n)) {
    rows = i+n*(seq_len(nodes)-1L)
    rule = mixture_weights_r(log_term[rows])
    if (is.nan(rule$value)) { out[i, ] = NaN; next }
    out[i, 1L] = rule$value
    keep = rule$weights > 0
    out[i, -1L] = drop(crossprod(rule$weights[keep], score[rows[keep], , drop = FALSE]))
  }
  out
}

mixture_moments_r = function(log_term, gradient, hessian = NULL, third = NULL, n, q = 0L) {
  d = ncol(gradient)
  out = list(value = rep(-Inf, n), gradient = matrix(0, n, d),
    hessian = if (!is.null(hessian)) array(0, c(n, q, d)) else NULL,
    third = if (!is.null(third)) array(0, c(n, q, q, d)) else NULL)
  nodes = if (n) length(log_term)/n else 0L
  for (i in seq_len(n)) {
    rows = i+n*(seq_len(nodes)-1L)
    rule = mixture_weights_r(log_term[rows])
    out$value[i] = rule$value
    if (is.nan(rule$value)) {
      out$gradient[i, ] = NaN
      if (!is.null(hessian)) out$hessian[i, , ] = NaN
      if (!is.null(third)) out$third[i, , , ] = NaN
      next
    }
    keep = rule$weights > 0
    if (!any(keep)) next
    rows = rows[keep]; w = rule$weights[keep]
    g = gradient[rows, , drop = FALSE]
    mean_g = drop(crossprod(w, g))
    out$gradient[i, ] = mean_g
    centered = sweep(g, 2L, mean_g)
    if (is.null(hessian)) next
    mean_h = matrix(0, q, d)
    for (j in seq_len(q)) for (k in seq_len(d)) {
      mean_h[j, k] = sum(w*hessian[rows, j, k])
      out$hessian[i, j, k] = mean_h[j, k] + sum(w*centered[, j]*centered[, k])
    }
    if (is.null(third)) next
    for (j in seq_len(q)) for (l in seq_len(q)) for (k in seq_len(d)) {
      out$third[i, j, l, k] = sum(w*third[rows, j, l, k]) +
        sum(w*(hessian[rows, j, l]-mean_h[j, l])*centered[, k]) +
        sum(w*(hessian[rows, j, k]-mean_h[j, k])*centered[, l]) +
        sum(w*(hessian[rows, l, k]-mean_h[l, k])*centered[, j]) +
        sum(w*centered[, j]*centered[, l]*centered[, k])
    }
  }
  out
}

affine_batch_r = function(base, location, transform, location_gradient = NULL,
                          transform_gradient = NULL) {
  groups = nrow(location); q = ncol(location); nodes = nrow(base)
  p = if (is.null(location_gradient)) 0L else dim(location_gradient)[3L]
  out = list(nodes = matrix(0, groups*nodes, q),
    gradient = if (!is.null(location_gradient)) array(0, c(groups*nodes, q, p)) else NULL)
  for (b in seq_len(nodes)) for (g in seq_len(groups)) {
    row = g+groups*(b-1L)
    out$nodes[row, ] = location[g, ] + matrix(transform[g, , ], q, q) %*% base[b, ]
    for (k in seq_len(p)) out$gradient[row, , k] = location_gradient[g, , k] +
      matrix(transform_gradient[g, , , k], q, q) %*% base[b, ]
  }
  out
}

row_contract_r = function(score, derivative) {
  n = nrow(score); q = ncol(score); p = dim(derivative)[3L]
  out = matrix(0, n, p)
  for (k in seq_len(p)) out[, k] = rowSums(score*matrix(derivative[, , k], n, q))
  out
}
