test_that('model declarations own ordered parameters and capabilities', {
  cases = list(
    list('probit', FALSE, 'panel', 'ProbitRE.R', 'sd_re_outcome'),
    list('poisson', FALSE, 'panel', 'PoissonRE.R', 'sd_re_outcome'),
    list('pln', FALSE, 'panel', 'PLNRE.R', c('sd_re_outcome', 'sd_error')),
    list('probit', TRUE, 'panel', 'ProbitRE_ProbitRE.R',
         c('sd_re_selection', 'sd_re_outcome', 'cor_re', 'cor_error')),
    list('gaussian', TRUE, 'panel', 'ProbitRE_LinearRE.R',
         c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error')),
    list('poisson', TRUE, 'panel', 'ProbitRE_PoissonRE.R',
         c('sd_re_selection', 'sd_re_outcome', 'cor_re')),
    list('pln', TRUE, 'panel', 'ProbitRE_PLNRE.R',
         c('sd_re_selection', 'sd_re_outcome', 'sd_error', 'cor_re', 'cor_error')),
    list('probit', TRUE, 'cross_sectional', 'Probit_Probit.R', 'cor_error'),
    list('gaussian', TRUE, 'cross_sectional', 'Probit_Linear.R', c('sd_error', 'cor_error')),
    list('pln', TRUE, 'cross_sectional', 'Probit_PLN.R', c('sd_error', 'cor_error')))
  for (case in cases) {
    model = resolve_model(case[[1]], case[[2]], case[[3]])
    expect_identical(model$file, case[[4]])
    expect_identical(model$parameters, case[[5]])
    expect_identical(model$capabilities$history, case[[3]] == 'panel')
    expect_true(is.function(model$panel) && is.function(model$predictions))
  }
  expect_error(resolve_model('probit', FALSE, 'cross_sectional'), 'require selection')
  expect_error(resolve_model('poisson', TRUE, 'cross_sectional'), 'requires a panel ID')
  expect_error(resolve_model(model_id = 'probit_poisson'), 'requires a panel ID')
  expect_identical(resolve_model('gaussian', TRUE)$positive_scales, 'sd_error')
  expect_identical(resolve_model('pln', TRUE, 'cross_sectional')$zero_correlations,
                   list(cor_error = 'sd_error'))
})

test_that('derived components retain unit counts and resolve their own model identity', {
  fixture = core_fixture('selection_pln')
  prepared = prepare_integration(fixture$prepared, core_control())
  expect_identical(prepared$model_id, 'probit_re_pln_re')
  expect_identical(prepared$integration_cache$selection$model_id, 'probit_re')
  expect_identical(prepared$integration_cache$outcome$model_id, 'pln_re')
  subset = subset_panels(prepared, 1L)
  expect_identical(subset$n_units, subset$n_panels)
  expect_null(subset$integration_cache)
})
