# Formula handling stays in R; the numerical evaluator receives aligned arrays.
model_design = function(formula, data, design = NULL, equation = "outcome") {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("The ", equation, " equation must be a two-sided formula", call. = FALSE)
  frame = if (is.null(design))
    stats::model.frame(formula, data, na.action = stats::na.pass, drop.unused.levels = TRUE) else
    stats::model.frame(design$terms, strip_contrasts(data), na.action = stats::na.pass,
                       drop.unused.levels = FALSE, xlev = design$xlevels)
  if (nrow(frame) != nrow(data))
    stop("Each equation must have one row per data observation", call. = FALSE)
  response = stats::model.response(frame)
  if (!is.numeric(response) && !is.logical(response))
    stop("Responses must be numeric or logical", call. = FALSE)
  if (is.matrix(response)) stop("Matrix responses are not supported", call. = FALSE)
  terms = attr(frame, "terms")
  matrix = stats::model.matrix(terms, frame, contrasts.arg = design$contrasts)
  if (!is.null(design) && !identical(colnames(matrix), design$columns))
    stop("Model columns do not match the fitted design", call. = FALSE)
  offset = stats::model.offset(frame)
  if (is.null(offset)) offset = rep(0, nrow(data))
  list(response = as.numeric(response), matrix = matrix, offset = offset, frame = frame,
       design = list(formula = formula, terms = terms,
                     xlevels = stats::.getXlevels(terms, frame),
                     contrasts = attr(matrix, "contrasts"), columns = colnames(matrix),
                     assign = attr(matrix, "assign")))
}

# model.frame(xlev =) re-levels factors and warns when that drops a contrasts
# attribute. The fitted design supplies the contrasts, so strip the attribute.
strip_contrasts = function(data) {
  for (j in seq_along(data))
    if (is.factor(data[[j]]) && !is.null(attr(data[[j]], "contrasts"))) attr(data[[j]], "contrasts") = NULL
  data
}

# A factor or character level present only in omitted rows would leave an
# all-zero design column. Re-encode the retained rows of the first-pass frame,
# so data-dependent terms (poly, ns, scale, I()) keep their evaluated values.
retain_levels = function(part, keep) {
  frame = part$frame
  categorical = which(vapply(frame, function(x) is.factor(x) || is.character(x), logical(1)))
  lost = vapply(categorical, function(j) {
    x = frame[[j]]
    used = if (is.factor(x)) levels(x) else unique(x[!is.na(x)])
    !all(used %in% x[keep])
  }, logical(1))
  if (!any(lost)) return(part)
  terms = attr(frame, "terms")
  retained = frame[keep, , drop = FALSE]
  for (j in categorical) if (is.factor(retained[[j]])) retained[[j]] = droplevels(retained[[j]])
  # A terms attribute makes model.matrix use the evaluated columns as they are.
  attr(retained, "terms") = terms
  matrix = stats::model.matrix(terms, retained)
  full = matrix(NA_real_, nrow(frame), ncol(matrix), dimnames = list(rownames(frame), colnames(matrix)))
  full[keep, ] = matrix
  part$matrix = full
  part$design$xlevels = stats::.getXlevels(terms, retained)
  part$design$contrasts = attr(matrix, "contrasts")
  part$design$columns = colnames(matrix)
  part$design$assign = attr(matrix, "assign")
  part
}

finite_rows = function(x) {
  if (!ncol(x)) return(rep(TRUE, nrow(x)))
  rowSums(!is.finite(x)) == 0L
}

# `.` would silently pull every other column (ID, responses) into an equation.
reject_dot = function(formula, equation) {
  if (inherits(formula, "formula") && "." %in% all.vars(formula))
    stop("`.` is not supported in the ", equation, " formula; list the regressors ",
         "(in the selection formula, `.outcome` stands for the outcome regressors)", call. = FALSE)
}

# `.outcome` in the selection formula stands for the outcome regressors: the
# outcome's term labels, without its offsets or its intercept setting. It may be
# combined by formula operators (`.outcome + z`, `.outcome:z`), not passed to a
# function such as log() or poly().
expand_outcome_terms = function(selection, outcome) {
  if (!inherits(selection, "formula") || length(selection) != 3L ||
      !".outcome" %in% all.vars(selection)) return(selection)
  if (".outcome" %in% all.vars(selection[[2L]]))
    stop("`.outcome` stands for regressors; it cannot be the selection response", call. = FALSE)
  labels = attr(stats::terms(outcome), "term.labels")
  if (!length(labels))
    stop("`.outcome` stands for the outcome regressors, but the outcome formula has none",
         call. = FALSE)
  regressors = str2lang(paste(labels, collapse = " + "))
  operators = c("+", "-", "*", ":", "/", "^", "%in%", "(")
  expand = function(e) {
    if (identical(e, quote(.outcome))) return(call("(", regressors))
    if (!is.call(e) || !".outcome" %in% all.vars(e)) return(e)
    if (!as.character(e[[1L]])[1L] %in% operators)
      stop("`.outcome` must be a term of the selection formula, not an argument of ",
           deparse(e[[1L]]), "()", call. = FALSE)
    for (i in seq_along(e)[-1L]) e[[i]] = expand(e[[i]])
    e
  }
  selection[[3L]] = expand(selection[[3L]])
  selection
}

# Rows selected by `subset =` as a logical vector over data rows: NA marks a
# missing row (na_action decides). NULL keeps every row.
subset_rows = function(value, data) {
  if (is.null(value)) return(NULL)
  n = nrow(data)
  if (is.logical(value)) {
    if (length(value) != n || !is.null(dim(value)))
      stop("A logical subset needs one value per data row", call. = FALSE)
    return(as.vector(value))
  }
  index = if (is.character(value)) match(value, rownames(data)) else if (is.numeric(value) &&
      !anyNA(value) && all(value == floor(value)) && all(abs(value) <= n) && all(value != 0) &&
      (all(value > 0) || all(value < 0))) seq_len(n)[value] else
    stop("subset must be logical, positive or negative row indices, or row names", call. = FALSE)
  if (anyNA(index)) stop("subset names rows that are not in data", call. = FALSE)
  if (anyDuplicated(index)) stop("subset must not repeat rows", call. = FALSE)
  seq_len(n) %in% index
}

prepare_data = function(formula, selection, data, id, family, na_action = "fail", design = NULL,
                        subset = NULL) {
  data = as.data.frame(data)
  if (!nrow(data)) stop("data must contain at least one row", call. = FALSE)
  family = family_name(family)
  has_selection = !is.null(selection)
  shorthand = NULL
  if (is.null(design)) {
    reject_dot(formula, "outcome")
    if (inherits(formula, "formula") && ".outcome" %in% all.vars(formula))
      stop("`.outcome` may appear only in the selection formula", call. = FALSE)
    if (has_selection) {
      reject_dot(selection, "selection")
      expanded = expand_outcome_terms(selection, formula)
      if (!identical(expanded, selection)) shorthand = selection
      selection = expanded
    }
  }
  cross_sectional = is.null(id)
  if (cross_sectional && !has_selection)
    stop("Cross-sectional estimation requires a selection equation; selection = NULL fits ",
         "only panel random-effects models (with an id column)", call. = FALSE)
  if (!cross_sectional) {
    if (!is.character(id) || length(id) != 1L || !id %in% names(data))
      stop("id must name one column in data", call. = FALSE)
    if (is.list(data[[id]]) || is.matrix(data[[id]]))
      stop("The id column must be an atomic vector", call. = FALSE)
  }
  structure = if (cross_sectional) "cross_sectional" else "panel"
  if (!has_selection && family == "gaussian")
    stop("The linear outcome model requires a selection equation; selection = NULL fits ",
         "only probit, poisson or pln random-effects models", call. = FALSE)
  model = resolve_model(family, has_selection, structure)
  ids_input = if (cross_sectional) seq_len(nrow(data)) else data[[id]]
  outcome = model_design(formula, data, design$outcome)
  select = if (has_selection) model_design(selection, data, design$selection, "selection") else NULL
  n = nrow(data)
  # Rows outside a subset are neither used nor checked; an NA subset entry is a
  # missing row.
  in_sample = if (is.null(subset)) rep(TRUE, n) else is.na(subset) | subset
  y = outcome$response
  d = if (has_selection) select$response else rep(1, n)
  if (any(in_sample & is.finite(d) & !d %in% c(0, 1)))
    stop("The selection response must be zero or one", call. = FALSE)
  selected = !is.na(d) & d == 1
  checked = selected & in_sample
  if (family == "probit" && any(is.finite(y[checked]) & !y[checked] %in% c(0, 1)))
    stop("Binary outcomes must be zero or one", call. = FALSE)
  if (family %in% c("poisson", "pln") &&
      any(is.finite(y[checked]) & (y[checked] < 0 | y[checked] != floor(y[checked]))))
    stop("Count outcomes must be nonnegative integers", call. = FALSE)

  valid = !is.na(ids_input) & is.finite(d)
  if (is.numeric(ids_input)) valid = valid & is.finite(ids_input)
  if (has_selection)
    valid = valid & finite_rows(select$matrix) & is.finite(select$offset)
  valid = valid & (!selected | (is.finite(y) & finite_rows(outcome$matrix) &
                                is.finite(outcome$offset)))
  if (!is.null(subset)) valid = valid & !is.na(subset)
  omitted = unname(which(!valid & in_sample))
  if (length(omitted) && na_action == "fail") {
    # Without a selection equation every outcome is used, so missing outcomes
    # often mean that the selection equation was left out.
    hint = if (!has_selection && any(!is.finite(y[omitted])))
      paste0("; outcomes may be missing only where a selection indicator is zero, so give ",
             "selection = indicator ~ regressors if they mark unselected rows") else
      " (na_action = \"omit\" or \"exclude\" drops such rows)"
    stop("Missing or nonfinite required data in row(s): ",
         paste(utils::head(omitted, 10L), collapse = ", "),
         if (length(omitted) > 10L) paste0(" and ", length(omitted) - 10L, " more"), hint, call. = FALSE)
  }
  keep = unname(which(valid & in_sample))
  if (!length(keep)) stop("No usable observations remain", call. = FALSE)
  if (is.null(design) && length(keep) < n) {
    outcome = retain_levels(outcome, keep)
    if (has_selection) select = retain_levels(select, keep)
  }
  ord = order(ids_input[keep], method = "radix")
  index = keep[ord]
  ids = ids_input[index]
  id_values = unique(ids)
  group = match(ids, id_values)
  selected = d[index] == 1
  x = outcome$matrix[index, , drop = FALSE]
  offset_outcome = outcome$offset[index]
  # Values never used by the likelihood must not propagate NA through arithmetic.
  x[!selected, ] = 0
  offset_outcome[!selected] = 0
  y = y[index]
  y[!selected] = 0
  omitted = if (length(omitted)) structure(omitted, names = rownames(data)[omitted],
                                          class = if (na_action == "exclude") "exclude" else "omit")
  prepared = list(model = if (cross_sectional) model$id else if (has_selection)
         paste0("selection_", family) else paste0(family, "_re"),
       model_id = model$id, structure = structure,
       family = family, has_selection = has_selection, id = id,
       x = x, w = if (has_selection) select$matrix[index, , drop = FALSE] else
         matrix(numeric(), length(index), 0L), y = y, d = d[index],
       offset_outcome = offset_outcome,
       offset_selection = if (has_selection) select$offset[index] else rep(0, length(index)),
       group = group, group_ptr = c(0L, cumsum(tabulate(group))), id_values = id_values,
       n_rows = length(index), n_units = length(id_values),
       n_panels = length(id_values), selected = selected,
       original_index = index, inverse_index = order(ord), keep = keep, na.action = omitted,
       design = list(outcome = outcome$design,
                     selection = if (has_selection) c(select$design,
                       if (!is.null(shorthand)) list(shorthand = shorthand)) else NULL))
  if (!is.null(subset)) prepared$subset = subset
  prepared
}

in_design_span = function(matrix, value) {
  if (all(value == 0)) return(TRUE)
  if (!ncol(matrix)) return(FALSE)
  unit = function(x) {
    scaled = x / max(abs(x))
    scaled / sqrt(sum(scaled^2))
  }
  matrix[] = vapply(seq_len(ncol(matrix)), function(j) unit(matrix[, j]),
                     numeric(nrow(matrix)))
  value = unit(value)
  coefficients = qr.coef(qr(matrix), value)
  scale = max(abs(coefficients), 1)
  scaled = coefficients / scale
  residual = sqrt(sum((value / scale - drop(matrix %*% scaled))^2))
  # Column normalization gives ||X||_F = sqrt(p); this is a backward-error test,
  # not a claim that every accepted design contains useful scale information.
  # Scale numerator and denominator together to avoid squaring large coefficients.
  backward = residual / (sqrt(ncol(matrix)) * sqrt(sum(scaled^2)) + 1 / scale)
  backward <= 64 * .Machine$double.eps * max(nrow(matrix), ncol(matrix) + 1L)
}

validate_singleton_selection = function(prepared, info, fixed) {
  if (identical(prepared$structure, "cross_sectional") ||
      !prepared$has_selection || prepared$n_panels != prepared$n_rows) return(invisible(NULL))
  free = function(name) name %in% info$name[!info$fixed]
  zero = function(name) name %in% names(fixed) && fixed[[name]] == 0
  index_can_scale = function(equation) {
    selection = equation == "selection"
    matrix = if (selection) prepared$w else prepared$x[prepared$selected, , drop = FALSE]
    known = if (selection) prepared$offset_selection else prepared$offset_outcome[prepared$selected]
    columns = if (ncol(matrix)) paste0(equation, ".", colnames(matrix)) else character()
    available = columns %in% info$name[!info$fixed]
    if (any(!available))
      known = known + drop(matrix[, !available, drop = FALSE] %*% fixed[columns[!available]])
    if (any(!is.finite(known)))
      stop("The fixed ", equation, " index is not finite; rescale the known predictors", call. = FALSE)
    in_design_span(matrix[, available, drop = FALSE], known)
  }
  selection_scale = free("sd_re_selection") && index_can_scale("selection")
  outcome_scale = if (prepared$family == "probit")
    free("sd_re_outcome") && index_can_scale("outcome") else
    prepared$family %in% c("gaussian", "pln") && free("sd_re_outcome") && free("sd_error")
  directions = c(sd_re_selection = selection_scale, sd_re_outcome = outcome_scale,
                  cor_re = free("cor_re"), cor_error = free("cor_error"))
  if (prepared$family %in% c("gaussian", "pln"))
    names(directions)[2L] = "sd_re_outcome/sd_error split"
  # Singleton models retain only one effective correlation. Two available
  # directions give a ridge; one does so only when its loadings are fixed zero.
  no_occasion = !"cor_error" %in% info$name || zero("cor_error")
  flat_selection = selection_scale && (zero("cor_re") || zero("sd_re_outcome")) &&
    (no_occasion || (prepared$family == "pln" && zero("sd_error")))
  flat_outcome = outcome_scale && (zero("cor_re") || zero("sd_re_selection")) && no_occasion
  flat_correlation = (free("cor_re") && (zero("sd_re_selection") || zero("sd_re_outcome"))) ||
    (free("cor_error") && prepared$family == "pln" && zero("sd_error"))
  if (sum(directions) > 1L || flat_selection || flat_outcome || flat_correlation)
    stop("A selection model with one retained row per panel has an unidentified singleton ",
         "likelihood ridge involving ", paste(names(directions)[directions], collapse = ", "),
         "; fix an appropriate scale or correlation, supply an independent known index ",
         "where applicable, supply repeated observations within panels, or, for independent ",
         "observations, use id = NULL (a cross-sectional model)", call. = FALSE)
}

validate_design_rank = function(prepared, info, fixed = numeric()) {
  if (prepared$model == "pln_re" && prepared$n_panels == prepared$n_rows &&
      all(!info$fixed[match(c("sd_re_outcome", "sd_error"), info$name)]))
    stop("A standalone PLN model with one retained row per panel cannot identify ",
         "both sd_re_outcome and sd_error; fix one scale or supply repeated ",
         "observations within panels", call. = FALSE)
  for (equation in c("selection", "outcome")) {
    if (equation == "selection" && !prepared$has_selection) next
    matrix = if (equation == "selection") prepared$w else
      prepared$x[prepared$selected, , drop = FALSE]
    columns = if (ncol(matrix)) paste0(equation, ".", colnames(matrix)) else character()
    free = !info$fixed[match(columns, info$name)]
    matrix = matrix[, free, drop = FALSE]
    if (ncol(matrix)) {
      decomposition = qr(matrix)
      if (!nrow(matrix) || decomposition$rank < ncol(matrix)) {
        # The QR pivot moves aliased columns behind the full-rank block.
        aliased = if (nrow(matrix))
          columns[free][decomposition$pivot[-seq_len(decomposition$rank)]] else character()
        stop("The free ", equation, " coefficients are not identified by the model matrix",
             if (length(aliased)) paste0(": ", paste(aliased, collapse = ", "),
               if (length(aliased) == 1L) " is" else " are",
               " aliased with other columns; drop them or fix them, e.g. fixed = c(`",
               aliased[1L], "` = 0)"), call. = FALSE)
      }
    }
    response = if (equation == "selection") prepared$d else prepared$y[prepared$selected]
    binary = equation == "selection" || prepared$family == "probit"
    separated = length(response) && if (binary) all(response == response[1L]) else
      prepared$family %in% c("poisson", "pln") && all(response == 0)
    # A free constant shift gives a strictly improving path to infinity.
    if (ncol(matrix) && separated &&
        max(abs(qr.resid(qr(matrix), rep(1, nrow(matrix))))) < 1e-10)
      stop("The ", equation, " response has no finite maximum-likelihood estimate ",
           "along a free constant predictor direction", call. = FALSE)
  }
  if (prepared$model == "probit_re" && prepared$n_panels == prepared$n_rows &&
      !info$fixed[match("sd_re_outcome", info$name)]) {
    columns = if (ncol(prepared$x)) paste0("outcome.", colnames(prepared$x)) else character()
    free = !info$fixed[match(columns, info$name)]
    known = prepared$offset_outcome
    if (any(!free))
      known = known + drop(prepared$x[, !free, drop = FALSE] %*% fixed[columns[!free]])
    if (any(!is.finite(known)))
      stop("The fixed probit index is not finite; rescale the known predictors",
           call. = FALSE)
    if (in_design_span(prepared$x[, free, drop = FALSE], known))
      stop("A standalone probit model with one row per panel cannot identify sd_re_outcome ",
           "when the known index is in the span of the free predictors at roundoff ",
           "tolerance; fix sd_re_outcome or supply an independent known index or repeated ",
           "observations within panels", call. = FALSE)
  }
  if (!any(prepared$selected) &&
      any(!info$fixed & (startsWith(info$name, "outcome.") |
                         info$name %in% c("sd_re_outcome", "sd_error", "cor_re", "cor_error"))))
    stop("No outcomes are observed; outcome parameters must be fixed", call. = FALSE)
  model = resolve_model(prepared$family, prepared$has_selection,
    prepared$structure, prepared$model_id)
  if (!is.null(model$validate_design)) model$validate_design(prepared, info)
  validate_singleton_selection(prepared, info, fixed)
}
