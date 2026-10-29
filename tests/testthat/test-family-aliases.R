test_that('family names and unique prefixes select the intended family', {
  accepted = list(
    gaussian = c('linear', 'continuous', 'gaussian', 'l', 'lin', 'cont', 'g', 'gauss'),
    probit = c('probit', 'binary', 'pr', 'prob', 'b', 'bin'),
    poisson = c('poisson', 'po', 'pois'),
    pln = c('pln', 'count', 'pl', 'cou'))
  for (family in names(accepted)) for (name in accepted[[family]])
    expect_identical(family_name(name), family, info = name)
  for (name in c('', 'p', 'c', 'co', 'count_p', 'count_pln', 'count_poisson',
    'count_full', 'normal', NA_character_))
    expect_error(family_name(name), 'is not available', info = name)
  expect_error(family_name(c('linear', 'probit')), 'family must be one of')
  expect_error(family_name(stats::poisson()), 'a family object')
})

test_that('aliases and abbreviations preserve all model structures and fitted quantities', {
  aliases = list(linear = c('continuous', 'gaussian'), probit = 'binary',
    poisson = 'pois', pln = 'count')
  data = data.frame(id = rep(1:3, each = 2), d = c(1, 1, 1, 0, 0, 1),
    x = c(-.4, .3, .7, -.2, .1, .5), z = c(.2, -.5, .1, .4, -.3, .6))
  outcome = y ~ x
  control = panel_control(nodes = 4L, observation_nodes = 6L)
  for (family in names(aliases)) for (structure in c('panel', 'cross_sectional', 'standalone')) {
    if (family == 'linear' && structure == 'standalone') next
    if (family == 'poisson' && structure == 'cross_sectional') next
    current = data
    current$y = switch(family, linear = c(.2, -.1, .8, .4, -.2, .5),
      probit = c(1, 0, 1, 0, 1, 0), c(1, 0, 2, 1, 0, 3))
    selection = if (structure != 'standalone') d ~ z
    id = if (structure != 'cross_sectional') 'id'
    fixed = c('outcome.(Intercept)' = .25, outcome.x = -.2)
    if (!is.null(selection)) {
      current$y[current$d == 0] = NA_real_
      fixed = c('selection.(Intercept)' = -.15, selection.z = .3, fixed)
    }
    if (!is.null(id)) {
      fixed = c(fixed, sd_re_outcome = .4)
      if (!is.null(selection)) fixed = c(fixed, sd_re_selection = .45, cor_re = -.35)
    }
    if (family %in% c('linear', 'pln')) fixed = c(fixed, sd_error = .5)
    if (!is.null(selection) && family != 'poisson') fixed = c(fixed, cor_error = .25)
    fit = function(name) panel_select(outcome, selection, current, id, name,
      fixed = fixed, control = control)
    reference = fit(family)
    for (alias in aliases[[family]]) {
      actual = fit(alias)
      expect_identical(actual$prepared, reference$prepared, info = paste(alias, structure))
      expect_identical(coef(actual), coef(reference))
      expect_identical(actual$loglik_by_unit, reference$loglik_by_unit)
      expect_identical(actual$score, reference$score)
      expect_identical(predict(actual, current), predict(reference, current))
      # Direct preparation and update() share the same normalization.
      expect_identical(prepare_data(y ~ x, selection, current, id, alias)$model_id,
        reference$model_id)
      restored = unserialize(serialize(actual, NULL))
      expect_identical(update(restored)$loglik_by_unit, reference$loglik_by_unit)
      expect_identical(update(reference, family = alias)$loglik_by_unit, reference$loglik_by_unit)
      expect_identical(update(reference, family = alias, evaluate = FALSE)$family, alias)
    }
  }
})

test_that('Poisson selection requires a panel ID on fitting and update', {
  data = data.frame(id = c(1, 1, 2, 2), d = c(1, 0, 1, 1), y = c(0, NA, 2, 1))
  for (family in c('poisson', 'pois', 'po'))
    expect_error(panel_select(y ~ 1, d ~ 1, data, NULL, family),
      'requires a panel ID.*use family = "count" or "pln"')
  expect_error(prepare_data(y ~ 1, d ~ 1, data, NULL, 'poisson'), 'requires a panel ID')
  fixed = c('selection.(Intercept)' = .2, 'outcome.(Intercept)' = .3)
  panel = panel_select(y ~ 1, d ~ 1, data, 'id', 'poisson',
    fixed = c(fixed, sd_re_selection = .4, sd_re_outcome = .5, cor_re = .2))
  expect_error(update(panel, id = NULL), 'requires a panel ID')
  cross = panel_select(y ~ 1, d ~ 1, data, NULL, 'count',
    fixed = c(fixed, sd_error = .5, cor_error = .2))
  expect_identical(cross$model_id, 'probit_pln')
  expect_error(update(cross, family = 'poisson', fixed = fixed), 'requires a panel ID')
})

test_that('aliases keep outcome validation and model restrictions', {
  data = data.frame(id = rep(1:3, each = 2), d = c(1, 1, 0, 1, 0, 1),
    y = c(.5, 1, NA, 2, NA, 0))
  for (family in c('poisson', 'count'))
    expect_error(panel_select(y ~ 1, d ~ 1, data, 'id', family),
      'Count outcomes must be nonnegative integers')
  expect_error(panel_select(y ~ 1, NULL, data, 'id', 'continuous'),
    'linear outcome model requires a selection equation')
  for (family in c('poisson', 'count'))
    expect_error(panel_select(y ~ 1, NULL, data, NULL, family),
      'requires a selection equation')
})
