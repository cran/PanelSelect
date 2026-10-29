# Numeric blocks use observations/panels on the first axis. Mixture and outer
# nodes are stacked node-major, with observation/panel varying fastest.
# Ordinary products retain IEEE arithmetic (including 0 * Inf = NaN).
group_sum = function(x, ptr) group_sum_cpp(x, ptr)

# Panel-fast rows, design columns plus a final constant one. Supplying normalized
# panel/node weights integrates local scores before expanding design columns.
group_weighted_design = function(x, weights, ptr, mixture_weights = NULL) {
  group_weighted_design_cpp(x, weights, ptr,
    if (is.null(mixture_weights)) matrix(numeric(), 0L, 0L) else mixture_weights,
    !is.null(mixture_weights))
}

logsum_score = function(log_term, score, n) logsum_score_cpp(log_term, score, n)

# Scalar log(sum(exp(x))) with a nonfinite maximum returned as is.
log_sum_exp = function(x) {
  peak = max(x)
  if (!is.finite(peak)) return(peak)
  peak + log(sum(exp(x-peak)))
}

logsum_weights = function(log_term, n) logsum_weights_cpp(log_term, n)

# R owns array axes; C++ sees ordinary matrices with flattened trailing axes.
numeric_block = function(x, shape) {
  if (!identical(dim(x), as.integer(shape))) stop('Array dimensions do not agree')
  matrix(x, shape[1L], prod(shape[-1L]))
}

# The C++ helper returns numerical moments, not model derivatives. Here the
# derivative of a log integral is assembled from its normalized node weights.
assemble_mixture_moments = function(moments, second = FALSE, third = FALSE) {
  list(value = moments$value, gradient = moments$mean_gradient,
    hessian = if (second) moments$mean_hessian + moments$cov_gradient else NULL,
    third = if (third) moments$mean_third + moments$cov_latent_hessian_score +
      moments$cov_mixed_hessian_score +
      aperm(moments$cov_mixed_hessian_score, c(1L, 3L, 2L, 4L)) +
      moments$central_third else NULL)
}

# Leading q gradient coordinates index H[, j, k] and T[, j, l, k].
mixture_moments = function(log_term, gradient, hessian = NULL, third = NULL, n, q = 0L) {
  m = nrow(gradient); d = ncol(gradient)
  if (length(q) != 1L || !is.finite(q) || q < 0 || q != floor(q))
    stop('Leading derivative dimension must be a nonnegative integer')
  if (q > d) stop('Leading derivative dimension exceeds the coordinate dimension')
  if (!is.null(third) && is.null(hessian)) stop('Third moments require second partials')
  if (is.null(hessian)) {
    out = logsum_score(log_term, gradient, n)
    return(list(value = out[, 1L], gradient = out[, -1L, drop = FALSE],
      hessian = NULL, third = NULL))
  }
  moments = mixture_moments_cpp(log_term, gradient,
    numeric_block(hessian, c(m, q, d)),
    if (is.null(third)) matrix(numeric(), m, 0L) else numeric_block(third, c(m, q, q, d)), n, q)
  for (field in c('mean_hessian', 'cov_gradient')) dim(moments[[field]]) = c(n, q, d)
  if (!is.null(third)) for (field in c('mean_third', 'cov_latent_hessian_score',
    'cov_mixed_hessian_score', 'central_third')) dim(moments[[field]]) = c(n, q, q, d)
  assemble_mixture_moments(moments, !is.null(hessian), !is.null(third))
}

# Apply all R-computed transforms with two BLAS products, preserving panel-fast
# rows within each node and (row, latent coordinate, parameter) derivatives.
affine_batch = function(base, location, transform, location_gradient = NULL,
                        transform_gradient = NULL) {
  groups = nrow(location); q = ncol(location); nodes = nrow(base)
  supplied = !is.null(location_gradient)
  if (supplied != !is.null(transform_gradient)) stop('Both affine derivative blocks are required')
  if (ncol(base) != q || !identical(dim(transform), as.integer(c(groups, q, q))))
    stop('Affine dimensions do not agree')
  index = rep(seq_len(groups), nodes)
  z = matrix(transform, groups*q, q) %*% t(base)
  z = matrix(aperm(array(z, c(groups, q, nodes)), c(1L, 3L, 2L)), groups*nodes, q)
  out = list(nodes = z + location[index, , drop = FALSE], gradient = NULL)
  if (supplied) {
    p = dim(location_gradient)[3L]
    if (!identical(dim(location_gradient), as.integer(c(groups, q, p))) ||
        !identical(dim(transform_gradient), as.integer(c(groups, q, q, p))))
      stop('Affine derivative dimensions do not agree')
    dz = matrix(aperm(transform_gradient, c(1L, 2L, 4L, 3L)), groups*q*p, q) %*% t(base)
    out$gradient = array(aperm(array(dz, c(groups, q, p, nodes)), c(1L, 4L, 2L, 3L)),
      c(groups*nodes, q, p)) + location_gradient[index, , , drop = FALSE]
  }
  out
}

row_contract = function(score, derivative) {
  p = dim(derivative)[3L]
  row_contract_cpp(score, numeric_block(derivative, c(nrow(score), ncol(score), p)), p)
}

