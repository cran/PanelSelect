# From observation contributions to panel derivatives.
#
# Read a model's panel function alongside add_rows: the model builds A and B,
# evaluates each row kind, and declares its natural first-score chains. Here,
# node_rows pairs observations with nodes and add_rows sums within each person.
# chain_score applies the declared chains only to the natural first scores;
# latent and higher mixed derivatives use the index maps from person_effects.
# Keep these two representations of the same mathematics consistent.
#
# make_panel_kernel binds one model point for the integration engines. It calls
# the model's panel function, except that standalone values and first scores
# use evaluate_standalone_nodes. It also adds the exact Linear outcome marginal
# and, when requested, the standard-normal log density of the nodes.
#
# Layout: z has one row per (person, node), person fastest within each node.
# q is the number of active latent coordinates; nq = q when latent derivatives
# are requested and 0 otherwise. p is the number of natural parameter columns,
# or 0 when parameters = FALSE. With N = nrow(z), the returned arrays are
# value[N], gradient[N, nq+p], hessian[N, nq, nq+p], third[N, nq, nq, nq+p].
# Hessian and third arrays are NULL unless their derivative order is requested.
# Latent columns come first; these are the mixed blocks needed for adaptive GH,
# not full parameter Hessians or third derivatives.

# 1. Compact local derivatives --------------------------------------------
# Each row is one observation/node entry. Of d local coordinates, the first
# r are predictor indices (a and/or b); remaining coordinates are natural
# error parameters. Arrays retain only the leading predictor axes:
# gradient[n, d], hessian[n, r, d], third[n, r, r, d]. A row function that
# already holds its first partials passes them as gradient (n x d, by column)
# rather than filling a zero matrix.
local_block = function(value, coordinates, r, order, gradient = NULL) {
  n = length(value); d = length(coordinates)
  if (!order) gradient = matrix(0, n, 0L)
  else if (is.null(gradient)) gradient = matrix(0, n, d)
  else dim(gradient) = c(n, d)
  list(value = value, gradient = gradient,
       hessian = if (order >= 2L) array(0, c(n, r, d)) else NULL,
       third = if (order >= 3L) array(0, c(n, r, r, d)) else NULL,
       coordinates = coordinates, r = r)
}

local_add = function(x, y) {
  for (name in c('value', 'gradient', 'hessian', 'third'))
    if (!is.null(x[[name]])) x[[name]] = x[[name]] + y[[name]]
  x
}

# 2. Person effects at the nodes ------------------------------------------

node_derivative_block = function(rows, dimension, nq, order) {
  list(value = numeric(rows), gradient = matrix(0, rows, dimension),
       hessian = if (order >= 2L) array(0, c(rows, nq, dimension)) else NULL,
       third = if (order >= 3L) array(0, c(rows, nq, nq, dimension)) else NULL)
}

# An effect location + loading*z (make_model_context) at the nodes, with its
# derivatives: latent columns when !compact, then the map's active parameter
# columns (compact keeps only those).
random_effect_values = function(map, z, panels, nq, p, order, compact = FALSE) {
  q = ncol(z); dimension = nq+p
  value = map$location[panels]
  if (q) value = value+rowSums(map$loading[panels, , drop = FALSE]*z)
  columns = if (p && order) map$columns else integer()
  gradient = matrix(0, nrow(z), if (compact) length(columns) else dimension)
  if (order) {
    if (nq && !compact) gradient[, seq_len(nq)] = map$loading[panels, , drop = FALSE]
    if (length(columns)) {
      derivative = map$location_gradient[panels, columns, drop = FALSE]
      for (j in seq_len(q)) derivative = derivative+
        matrix(map$loading_gradient[panels, j, columns], nrow(z), length(columns))*z[, j]
      gradient[, if (compact) seq_along(columns) else nq+columns] = derivative
    }
  }
  list(value = value, gradient = gradient, columns = columns)
}

# Derivatives of the index scale*effect for the latent and mixed blocks: the
# latent loading, the parameter columns (including the effect's own SD) and the
# latent x parameter terms.
map_index_derivatives = function(effect, map, scale, scale_index, panels, q, nq, p, order) {
  if (!order) return(NULL)
  if (order == 1L)
    return(list(loading = if (nq) map$loading[panels, , drop = FALSE]*scale else NULL))
  groups = length(panels); parameter_columns = nq+seq_len(p)
  if (p) {
    extra = effect$gradient[, parameter_columns, drop = FALSE]*scale
    if (!is.na(scale_index)) extra[, scale_index] = extra[, scale_index]+effect$value
    mixed = map$loading_gradient[panels, , , drop = FALSE]*scale
    dim(mixed) = c(groups, q, p)
    if (!is.na(scale_index)) mixed[, , scale_index] =
      matrix(mixed[, , scale_index], groups, q)+map$loading[panels, , drop = FALSE]
  } else {
    extra = matrix(0, groups, p)
    mixed = array(0, c(groups, q, p))
  }
  loading = map$loading[panels, , drop = FALSE]*scale
  list(extra = extra, columns = seq_len(p), loading = loading, mixed = mixed)
}

# Person effects at the nodes z: values u and v, their natural-parameter
# derivatives du and dv (v_rho = dv/drho), the index maps of A and B for latent
# columns, and the columns du and dv cover.
person_effects = function(m, z, out, effects = c('u', 'v')) {
  compact = out$order < 2L; nq = out$nq; p = out$p
  u = if ('u' %in% effects) random_effect_values(m$u, z, out$panels, nq, p, out$order, compact)
  v = if ('v' %in% effects) random_effect_values(m$v, z, out$panels, nq, p, out$order, compact)
  e = list(u = u$value, v = v$value, index = list(
    b = if (!is.null(v)) map_index_derivatives(v, m$v, m$lambda, m$lambda_index,
      out$panels, m$q, nq, p, out$order),
    a = if (!is.null(u)) map_index_derivatives(u, m$u, m$delta,
      match('sd_re_selection', out$names), out$panels, m$q, nq, p, out$order)),
    columns = list(u = nq+u$columns, v = nq+v$columns))
  if (p && out$order) {
    effect_gradient = function(effect) {
      if (is.null(effect) || !length(effect$columns)) return(NULL)
      if (compact) effect$gradient else effect$gradient[, nq+effect$columns, drop = FALSE]
    }
    e$du = effect_gradient(u); e$dv = effect_gradient(v)
    rho_column = match('cor_re', out$names)
    e$v_rho = if (!is.null(v) && rho_column %in% v$columns)
      e$dv[, match(rho_column, v$columns)] else numeric(out$groups)
  }
  e
}

# Rows of one row block at every node of the batch: the observation row, the
# (person, node) row of z and the response.
node_rows = function(out, block) {
  count = length(block$rows)
  list(row = rep(block$rows, out$blocks),
       unit = rep(block$group, out$blocks)+rep((seq_len(out$blocks)-1L)*out$n, each = count),
       y = rep(block$y, out$blocks))
}

# 3. Sums within persons and the chain rule -------------------------------

# Add one row block's log contributions and natural first scores, preserving
# the declared term order. A, B, error and tau describe the chains for local
# coordinates a, b, error and tau. For S = sum_t l_coordinate, a vector term
# contributes term*S, a matrix under role u or v contributes its map columns
# times S, and 1 contributes S directly. Character terms ('w' or 'x') mark a
# design product: the local coordinate chooses block$w for a or block$x for b;
# the string itself is not a design selector. Role names choose score columns.
# Higher mixed derivatives below use e$index and structural_columns instead
# of interpreting these first-score chains. Empty blocks return before the
# local calculation is forced, so missing outcomes are never evaluated.
add_rows = function(out, e, local, block, A = NULL, B = NULL, error = NULL, tau = NULL) {
  count = length(block$rows)
  if (!count) return(out)
  chain = list(a = A, b = B, error = error, tau = tau)
  groups = out$groups; blocks = out$blocks; nq = out$nq; p = out$p
  sum_columns = function(x) matrix(group_sum(matrix(x, count), block$ptr), groups)
  out$value = out$value+drop(sum_columns(local$value))
  if (!out$order) return(out)
  d = length(local$coordinates); r = local$r
  design = function(name) if (name == 'a') block$w else block$x
  # Single-index scores need one weighted design reduction. Its final column
  # is the score sum, so no separate grouping or full-width pullback is needed.
  if (out$order == 1L && p && d == 1L) {
    name = local$coordinates[1L]
    product = group_weighted_design(design(name), matrix(local$gradient, count, blocks), block$ptr)
    first = product[, ncol(product)]
    out$gradient = chain_score(out$gradient, chain[name], list(first), list(product), out, e)
    for (j in seq_len(nq)) out$gradient[, j] =
      out$gradient[, j]+first*e$index[[name]]$loading[, j]
    return(out)
  }
  first = sum_columns(local$gradient)
  if (p) {
    products = vector('list', d)
    for (i in seq_len(d)) {
      name = local$coordinates[i]
      if (name %in% c('a', 'b') && length(if (name == 'a') out$iw else out$ix))
        products[[i]] = group_weighted_design(design(name),
          matrix(local$gradient[, i], count, blocks), block$ptr)
    }
    out$gradient = chain_score(out$gradient, chain[local$coordinates],
      lapply(seq_len(d), function(i) first[, i]), products, out, e)
  }
  loading = lapply(local$coordinates[seq_len(r)], function(name) e$index[[name]]$loading)
  for (j in seq_len(nq)) for (i in seq_len(r))
    out$gradient[, j] = out$gradient[, j]+first[, i]*loading[[i]][, j]
  if (out$order < 2L || !nq) return(out)
  add_row_curvature(out, e, local, block, first, loading)
}

# Adds each local coordinate's declared terms to natural first-score columns.
# chains, sums and products align by coordinate: sums[[i]] is the partial
# summed within each person and products[[i]] its design product (or NULL).
# One call takes every coordinate of a row block, so the gradient matrix is
# copied once rather than once per coordinate.
chain_score = function(gradient, chains, sums, products, out, e) {
  for (i in seq_along(chains)) for (role in names(chains[[i]])) {
    term = chains[[i]][[role]]
    columns = if (role %in% c('u', 'v')) e$columns[[role]] else out$columns[[role]]
    if (!length(term) || !length(columns)) next
    gradient[, columns] = gradient[, columns, drop = FALSE]+
      if (is.character(term)) products[[i]][, seq_along(columns), drop = FALSE] else
      if (identical(term, 1)) sums[[i]] else term*sums[[i]]
  }
  gradient
}

# Latent Hessian and third derivatives with their mixed parameter blocks, for
# adaptive GH. The chain through A = a + delta*u(z) and B = b + lambda*v(z)
# comes from e$index; structural_columns maps error and tau to their natural
# parameters. This path does not consume the model's declared first-score chain.
add_row_curvature = function(out, e, local, block, first, loading) {
  count = length(block$rows); groups = out$groups; blocks = out$blocks
  nq = out$nq; p = out$p; d = length(local$coordinates); r = local$r
  parameter_columns = nq+seq_len(p)
  index_names = local$coordinates[seq_len(r)]
  sum_columns = function(x) matrix(group_sum(matrix(x, count), block$ptr), groups)
  pullback = function(derivative, sums)
    pullback_unit_derivatives(derivative, sums, local$coordinates, block,
      e$index, out, p, groups, blocks)
  second = array(sum_columns(local$hessian), c(groups, r, d))
  second_pullback = if (p) lapply(seq_len(r), function(i)
    pullback(matrix(local$hessian[, i, ], count*blocks, d), matrix(second[, i, ], groups, d)))
  for (j in seq_len(nq)) {
    if (p) {
      mixed = matrix(0, groups, p)
      for (i in seq_len(r)) {
        mixed = mixed+second_pullback[[i]]*loading[[i]][, j]+
          matrix(e$index[[index_names[i]]]$mixed[, j, ], groups, p)*first[, i]
      }
      out$hessian[, j, parameter_columns] =
        matrix(out$hessian[, j, parameter_columns], groups, p)+mixed
    }
    for (k in seq_len(nq)) for (i in seq_len(r)) for (l in seq_len(r))
      out$hessian[, j, k] = out$hessian[, j, k]+second[, i, l]*loading[[i]][, j]*loading[[l]][, k]
  }
  if (out$order < 3L) return(out)
  third = array(sum_columns(local$third), c(groups, r, r, d))
  third_pullback = if (p) lapply(seq_len(r*r), function(index) {
    i = (index-1L)%%r+1L; l = (index-1L)%/%r+1L
    pullback(matrix(local$third[, i, l, ], count*blocks, d), matrix(third[, i, l, ], groups, d))
  })
  for (j in seq_len(nq)) for (k in seq_len(nq)) {
    if (p) {
      mixed = matrix(0, groups, p)
      for (i in seq_len(r)) for (l in seq_len(r)) {
        mixed = mixed+third_pullback[[i+(l-1L)*r]]*(loading[[i]][, j]*loading[[l]][, k])+
          (matrix(e$index[[index_names[i]]]$mixed[, j, ], groups, p)*loading[[l]][, k]+
           matrix(e$index[[index_names[l]]]$mixed[, k, ], groups, p)*loading[[i]][, j])*second[, i, l]
      }
      out$third[, j, k, parameter_columns] =
        matrix(out$third[, j, k, parameter_columns], groups, p)+mixed
    }
    for (m in seq_len(nq)) for (i in seq_len(r)) for (l in seq_len(r)) for (s in seq_len(r))
      out$third[, j, k, m] = out$third[, j, k, m]+
        third[, i, l, s]*loading[[i]][, j]*loading[[l]][, k]*loading[[s]][, m]
  }
  out
}

# One slice of local higher derivatives expanded to natural parameters: design
# products for the indices a and b, the index maps' parameter columns, and the
# structural columns (sd_error for the error coordinate, cor_error for tau).
pullback_unit_derivatives = function(derivative, sums, coordinates, block,
                                      maps, context, p, groups, blocks) {
  d = length(coordinates); count = length(block$rows)
  ans = matrix(0, groups, p)
  for (i in seq_len(d)) {
    name = coordinates[i]
    if (name %in% c('a', 'b')) {
      map = maps[[name]]
      design = if (name == 'a') block$w else block$x
      columns = if (name == 'a') context$iw else context$ix
      if (length(columns)) {
        product = group_weighted_design(design,
          matrix(derivative[, i], count, blocks), block$ptr)
          ans[, columns] = ans[, columns, drop = FALSE]+product[, seq_along(columns), drop = FALSE]
      }
      if (length(map$columns)) ans[, map$columns] =
        ans[, map$columns, drop = FALSE]+map$extra*sums[, i]
    } else {
      column = unname(context$structural_columns[name])
      if (!is.na(column)) ans[, column] = ans[, column]+sums[, i]
    }
  }
  ans
}

# 4. Standalone first-order fast path -------------------------------------
# Value and first score of a standalone model (B = b + lambda*z) at order <= 1,
# with one design reduction per batch: the model's outcome rows, then
# g_beta = sum_t x_t l_B, g_lambda = z sum_t l_B and g_gamma = sum_t l_gamma.
# This avoids repeated panel setup in staged GL's per-panel evaluations while
# retaining the model's outcome rows and the shared first-score assembly.
evaluate_standalone_nodes = function(z, prior, prepared, context, order, p, nq) {
  q = context$q; n = prepared$n_panels
  block = context$data$active; count = length(block$rows)
  b = context$b[block$rows]; y = block$y
  factorial = context$data$log_factorial[block$rows]
  groups = nrow(z); blocks = groups/n
  out = node_derivative_block(groups, nq+p, nq, order)
  if (count) {
    panel = rep(block$group, blocks)+rep((seq_len(blocks)-1L)*n, each = count)
    B = rep(b, blocks)
    if (q) B = B+context$lambda*z[panel, 1L]
    local = context$definition$outcome_rows(B, rep(y, blocks), rep(factorial, blocks),
      context, order)
    value = local$value
    # array() keeps the leading count*blocks values: column 1 (the B partial)
    # of the column-major gradient, as rows x node blocks, in one copy.
    if (order) first = array(local$gradient, c(count, blocks))
    sum_columns = function(x) matrix(group_sum(matrix(x, count), block$ptr), groups)
    out$value = drop(sum_columns(value))
    if (order) {
      if (p) {
        product = group_weighted_design(block$x, first, block$ptr)
        total = product[, ncol(product)]
        chains = list(list(beta = 'x', lambda = if (q) z[, 1L] else numeric(groups)))
        sums = list(total)
        if ('error' %in% local$coordinates) {
          chains[[2L]] = list(gamma = 1); sums[[2L]] = drop(sum_columns(local$gradient[, 2L]))
        }
        out$gradient = chain_score(out$gradient, chains, sums, list(product, NULL),
          list(columns = context$score_columns), NULL)
      } else total = drop(sum_columns(first))
      if (nq) out$gradient[, 1L] = context$lambda*total
    }
  }
  if (prior && q) {
    out$value = out$value-log(2*pi)/2-z[, 1L]^2/2
    if (order && nq) out$gradient[, 1L] = out$gradient[, 1L]-z[, 1L]
  }
  out
}

# 5. Bind one parameter point ----------------------------------------------
# kernel$evaluate(z, prior) evaluates the model at z, using the standalone fast
# path for order <= 1 and the panel function otherwise. Shared contributions
# are the exact Linear marginal and optional normal prior on z; the prior is
# node-dependent. kernel$effects(z) returns u, v and their derivatives for
# history predictions, using the same model point and column layout.
make_panel_kernel = function(theta, prepared, control, fixed,
                             order = 1L, parameters = TRUE, latent = FALSE, context = NULL) {
  if (is.null(context)) context = make_model_context(theta, prepared, control, fixed)
  q = context$q; nq = if (latent) q else 0L
  p = if (parameters) length(theta) else 0L
  # Mode searches need index partials only; final geometry keeps all columns.
  context$parameter_derivatives = parameters
  if (p && nq) context$score_columns = lapply(context$score_columns, function(x) nq+x)
  panel = context$definition$panel
  evaluate = function(z, prior = FALSE) {
    recorder = context$staged_recorder %||% context$inner$recorder
    if (!is.null(recorder)) recorder$outer_radius = max(recorder$outer_radius, abs(z), 0)
    # Standalone values and first scores take the fast path of section 4.
    if (!prepared$has_selection && order <= 1L)
      return(evaluate_standalone_nodes(z, prior, prepared, context, order, p, nq))
    # The output arrays and the layout the shared helpers need to fill them.
    n = prepared$n_panels; groups = nrow(z)
    out = c(node_derivative_block(groups, nq+p, nq, order),
      list(order = order, nq = nq, p = p, n = n, groups = groups, blocks = groups/n,
        panels = rep(seq_len(n), groups/n), names = names(theta),
        columns = context$score_columns, iw = context$iw, ix = context$ix,
        structural_columns = context$structural_columns))
    out = panel(z, context, out)
    # The exact Linear outcome marginal is constant across a person's nodes.
    # The optional standard-normal log density below depends on z and adds
    # latent derivatives, but has no natural parameter derivatives at fixed z.
    if (!is.null(context$marginal)) {
      out$value = out$value+context$marginal$value[out$panels]
      if (order && p) out$gradient[, nq+seq_len(p)] =
        out$gradient[, nq+seq_len(p), drop = FALSE]+context$marginal$gradient[out$panels, , drop = FALSE]
    }
    if (prior && q) {
      out$value = out$value-q/2*log(2*pi)-rowSums(z*z)/2
      if (order && nq) out$gradient[, seq_len(nq)] = out$gradient[, seq_len(nq), drop = FALSE]-z
      if (order >= 2L) for (j in seq_len(nq)) out$hessian[, j, j] = out$hessian[, j, j]-1
    }
    out[c('value', 'gradient', 'hessian', 'third')]
  }
  effects = function(z) {
    panels = rep_len(seq_len(prepared$n_panels), nrow(z))
    list(u = random_effect_values(context$u, z, panels, nq, p, order),
         v = random_effect_values(context$v, z, panels, nq, p, order))
  }
  list(evaluate = evaluate, effects = effects, q = q,
       parameter_columns = nq+seq_len(p), latent_columns = seq_len(nq))
}
