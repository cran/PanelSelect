test_that('native contractions reject negative dimensions before reducing', {
  expect_error(row_contract_cpp(matrix(1, 2L, 1L), matrix(1, 2L, 1L), -1L), 'dimensions')
  expect_error(mixture_moments_cpp(c(0, 0), matrix(0, 2L, 1L), matrix(0, 2L, 0L),
    matrix(0, 2L, 0L), 1, -1L), 'dimensions')
  # Native weights return plain vectors to R.
  expect_null(dim(logsum_weights_cpp(c(0, log(2)), 1)$value))
})

test_that('cross-sectional rows share the panel row likelihood functions', {
  local_mocked_bindings(
    probit_re_probit_re_rows = function(...) stop('binary rows reached'),
    linear_selection_log_contribution = function(...) stop('Linear rows reached'),
    probit_re_pln_re_rows = function(...) stop('PLN rows reached'))
  rows = c(probit = 'binary rows reached', gaussian = 'Linear rows reached',
    pln = 'PLN rows reached')
  for (family in names(rows)) {
    fixture = core_fixture(paste0('selection_', family))
    fixture$prepared$structure = 'cross_sectional'
    fixture$prepared$model_id = c(probit = 'probit_probit', gaussian = 'probit_linear',
      pln = 'probit_pln')[[family]]
    theta = fixture$theta[!names(fixture$theta) %in% c('sd_re_selection', 'sd_re_outcome', 'cor_re')]
    expect_error(evaluate_model(theta, fixture$prepared, core_control(inner = 5L)),
      rows[[family]], fixed = TRUE)
  }
})

test_that('missing structure metadata is never read as a panel model', {
  fixture = core_fixture('selection_probit')
  stripped = fixture$prepared; stripped$model_id = stripped$structure = NULL
  expect_error(parameter_layout(stripped), 'missing metadata')
  expect_error(resolve_fixed(numeric(), parameter_layout(fixture$prepared), stripped),
    'missing metadata')
  # A recorded model identity fixes the structure of an older fitted object.
  expect_identical(model_structure(list(model_id = 'probit_pln', prepared = list())),
    'cross_sectional')
  expect_identical(model_structure(list(prepared = list(model_id = 'probit_re_pln_re'))), 'panel')
  expect_error(model_structure(list(prepared = list(family = 'probit'))), 'missing metadata')
})
