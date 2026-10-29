# update(), refit records and history rows: formula and selection updates,
# shorthand re-expansion, merged starts, lightweight call records, refits
# through update(), and history that follows the fit's missing-data rule.

update_fixture = function(seed = 1101) {
  local_test_seed(seed)
  panels = 20L; waves = 4L; n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n), z = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + stats::rnorm(n) > 0)
  data$y = stats::rpois(n, exp(.2 + .4*data$x + .6*v[data$id]))
  data$y[data$d == 0] = NA
  data
}

# Fits below start from pooled regressions (init = 'pooled').
update_control = function(...) panel_control(nodes = 4L, reltol = 1e-12, ...)

test_that('outcome and selection updates equal direct refits from the same warm starts', {
  data = update_fixture()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = update_control(), init = 'pooled')
  start = coef(fit)
  direct = function(outcome, selection, keep) panel_select(outcome,
    selection, data, 'id', 'poisson', start = start[intersect(names(start), keep)], control = update_control(), init = 'pooled')
  wider = update(fit, outcome = . ~ . + z)
  reference = direct(y ~ x + z, d ~ w, names(start))
  expect_identical(coef(wider), coef(reference))
  selection = update(fit, selection = . ~ . + x)
  expect_identical(coef(selection), coef(direct(y ~ x, d ~ w + x, names(start))))
  expect_equal(formula(selection, 'selection'), d ~ w + x, ignore_attr = TRUE)
  expect_error(update(panel_select(y ~ x, NULL, data[data$d == 1, ], 'id', 'poisson',
    control = update_control(), init = 'pooled'), selection = . ~ . + x), 'no selection equation to update')
  expect_error(update(fit, foo = 1), 'outcome, selection, data, id, family, start, fixed, control, subset, na_action, init, vcov')
})

test_that('a .outcome selection follows the updated outcome', {
  data = update_fixture()
  fit = panel_select(y ~ x + z, d ~ .outcome + w, data, 'id', 'poisson',
    control = update_control(), init = 'pooled')
  smaller = update(fit, . ~ . - z)
  expect_equal(formula(smaller, 'selection'), d ~ (x) + w, ignore_attr = TRUE)
  expect_equal(smaller$prepared$design$selection$shorthand, d ~ .outcome + w,
    ignore_attr = TRUE)
  expect_false(any(grepl('\\.z$', names(coef(smaller)))))
})

test_that('supplied starts merge with warm starts and NULL requests automatic starts', {
  data = update_fixture()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = update_control(), init = 'pooled')
  recipe = update(fit, start = c(outcome.x = .05), evaluate = FALSE)
  expected = coef(fit); expected[['outcome.x']] = .05
  expect_identical(recipe$start, expected)
  expect_null(update(fit, start = NULL, evaluate = FALSE)$start)
  expect_identical(names(as.list(recipe))[2L], 'outcome')
})

test_that('evaluated updates and refits store lightweight call records', {
  data = update_fixture()
  replacement = transform(data, x = x + .1)
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = update_control(), init = 'pooled')
  updated = update(fit)
  expect_false(is.data.frame(updated$call$data))
  expect_identical(updated$call$data, quote(data))
  expect_identical(update(fit, data = replacement)$call$data, quote(replacement))
  expect_lte(length(serialize(updated, NULL)), 1.05*length(serialize(fit, NULL)))
  # A record from an older version: formula named, data embedded, no explicit id.
  old = fit
  old$call = as.call(list(quote(panel_select), formula = y ~ x, data = data, id = 'id',
    family = 'poisson', selection = d ~ w))
  record = update(old)$call
  expect_identical(names(record)[2L], 'outcome')
  expect_identical(record$data, as.name('<retained data>'))
  cross = data; cross$id = NULL
  single = panel_select(y ~ x, d ~ w, cross, NULL, 'count',
    fixed = c(sd_error = 0, cor_error = 0), control = update_control(), init = 'pooled')
  single$call$id = NULL
  expect_true('id' %in% names(update(single)$call))
  expect_null(update(single)$call$id)
  refit = check_quadrature(fit, nodes = 6L, refit = TRUE)$refit
  expect_identical(refit$call$data, quote(data))
  expect_identical(refit$call$control$nodes, 6L)
})

test_that('quadrature refits through update() reproduce a direct warm-started refit', {
  data = update_fixture()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = update_control(), init = 'pooled')
  refined = fit$control; refined$nodes = 6L; refined = validate_panel_control(refined)
  refit = check_quadrature(fit, nodes = 6L, refit = TRUE)$refit
  direct = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', start = coef(fit),
    control = refined)
  expect_identical(coef(refit), coef(direct))
})

test_that('history rows follow the fit missing-data rule', {
  data = update_fixture()
  gappy = data; gappy$x[c(5L, 22L)] = NA
  fit = panel_select(y ~ x, d ~ w, gappy, 'id', 'poisson',
    control = update_control(), na_action = 'omit', init = 'pooled')
  targets = data[data$id %in% c(2, 6), ]
  expect_equal(predict(fit, targets, conditioning = 'history', history = gappy),
    predict(fit, targets, conditioning = 'history', history = gappy[!is.na(gappy$x), ]))
  expect_equal(random_effects(fit, gappy), random_effects(fit, gappy[!is.na(gappy$x), ]))
  strict = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = update_control(), init = 'pooled')
  expect_error(predict(strict, targets, conditioning = 'history', history = gappy),
    'Missing or nonfinite')
})
