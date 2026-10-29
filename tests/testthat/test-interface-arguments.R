# The fitting interface: outcome and selection come first, and data, id
# and family are required; the fit settings na_action, init and vcov are
# arguments of panel_select(), not of panel_control().

arguments_data = function(seed = 1201, panels = 20L, waves = 4L) {
  local_test_seed(seed)
  n = panels*waves
  data = data.frame(id = rep(seq_len(panels), each = waves), x = stats::rnorm(n),
    w = stats::rnorm(n))
  u = stats::rnorm(panels); v = .4*u + sqrt(1 - .4^2)*stats::rnorm(panels)
  data$d = as.integer(.3 + .6*data$w + .8*u[data$id] + stats::rnorm(n) > 0)
  data$y = stats::rpois(n, exp(.2 + .4*data$x + .6*v[data$id]))
  data$y[data$d == 0] = NA
  data$count = stats::rpois(n, exp(.1 + .3*data$x + .9*u[data$id]))
  data$b = as.integer(.2 + .5*data$x + .7*u[data$id] + stats::rnorm(n) > 0)
  data
}

arguments_control = function(...) panel_control(nodes = 4L, ...)

# Every component except the elapsed times.
fit_contents = function(fit) {
  fit$time = NULL; fit$optimizer$time = NULL
  unclass(fit)
}

test_that('positional arguments in the new order equal the fully named call', {
  data = arguments_data()
  positional = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = arguments_control(),
    init = 'pooled')
  named = panel_select(outcome = y ~ x, selection = d ~ w, data = data, id = 'id',
    family = 'poisson', control = arguments_control(), init = 'pooled')
  expect_identical(fit_contents(positional), fit_contents(named))
  standalone = panel_select(count ~ x, NULL, data, 'id', 'poisson', control = arguments_control())
  expect_identical(fit_contents(standalone), fit_contents(panel_select(outcome = count ~ x,
    selection = NULL, data = data, id = 'id', family = 'poisson', control = arguments_control())))
  expect_false(standalone$prepared$has_selection)
})

test_that('missing required arguments and the earlier order stop with the fix', {
  data = arguments_data()
  expect_error(panel_select(), 'outcome is required')
  expect_error(panel_select(count ~ x), 'selection = NULL for a standalone probit, Poisson or PLN')
  expect_error(panel_select(count ~ x, NULL), 'data is required')
  expect_error(panel_select(count ~ x, NULL, data), 'id = NULL for independent cross-sectional')
  expect_error(panel_select(count ~ x, NULL, data, 'id'),
    'family is required: "linear" (aliases "continuous" and "gaussian"), "probit" (alias "binary")',
    fixed = TRUE)
  # The earlier positional order put data second; naming data = shifts id and family left.
  expect_error(panel_select(count ~ x, data, 'id', 'poisson'),
    'takes \\(outcome, selection, data, id, family\\)')
  expect_error(panel_select(count ~ x, data = data, 'id', 'poisson'),
    'takes \\(outcome, selection, data, id, family\\).*received an object of class character')
  # NULL or a family object never selects a family silently.
  expect_error(panel_select(y ~ x, d ~ w, data, 'id', NULL), 'family is required')
  expect_error(panel_select(count ~ x, NULL, data, 'id', poisson()), 'a family object such as binomial')
  expect_error(panel_select(y ~ x, NULL, data, 'id', 'linear'),
    'requires a selection equation; selection = NULL fits only probit, poisson or pln')
  expect_error(panel_select(count ~ x, NULL, data, NULL, 'poisson'),
    'requires a selection equation; selection = NULL fits only panel random-effects models')
})

test_that('a forgotten selection equation is named when outcomes are missing', {
  data = arguments_data()
  expect_error(panel_select(y ~ x, NULL, data, 'id', 'poisson'),
    'outcomes may be missing only where a selection indicator is zero')
  # Missing regressors alone do not suggest a selection equation.
  gappy = data; gappy$x[3L] = NA
  message = tryCatch(panel_select(count ~ x, NULL, gappy, 'id', 'poisson'),
    error = conditionMessage)
  expect_match(message, 'row\\(s\\): 3 \\(na_action = "omit" or "exclude" drops such rows\\)$')
})

test_that('standalone models with every parameter fixed reproduce the frozen likelihoods', {
  for (model in c('probit_re', 'poisson_re', 'pln_re')) {
    reference = core_references$fixed[[paste(model, 1)]]
    family = sub('_re$', '', model)
    data = data.frame(id = 1L, x = reference$panel$x[, 'x'], y = reference$panel$y)
    fit = panel_select(y ~ x, NULL, data, 'id', family, fixed = reference$theta,
      control = panel_control(quadrature = 'gh', nodes = reference$outer,
        observation_nodes = reference$inner))
    expect_equal(fit$loglik, reference$loglik, tolerance = 2e-10, info = model)
    expect_identical(fit$optimizer$status, 'all_fixed')
  }
})

test_that('fit settings act as arguments of panel_select', {
  data = arguments_data()
  default = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = arguments_control(),
    init = 'pooled')
  observed = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = arguments_control(),
    init = 'pooled', vcov = 'observed')
  # add_information() recomputes the information at the reported estimates, a
  # second public path to the same matrix up to coordinate round trips.
  expect_equal(vcov(observed), vcov(add_information(default), type = 'observed'), tolerance = 1e-10)
  expect_identical(coef(observed), coef(default))
  expect_identical(observed$control$vcov, 'observed')
  gappy = data; gappy$x[c(2L, 30L)] = NA
  omitted = panel_select(y ~ x, d ~ w, gappy, 'id', 'poisson', control = arguments_control(),
    init = 'pooled', na_action = 'omit')
  complete = panel_select(y ~ x, d ~ w, gappy[!is.na(gappy$x), ], 'id', 'poisson',
    control = arguments_control(), init = 'pooled')
  expect_identical(coef(omitted), coef(complete))
  expect_identical(omitted$loglik, complete$loglik)
  # init = "zero" starts every optimizer coordinate at zero: coefficients 0, SDs 1
  # and correlations 0 on the natural scale.
  zero = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = arguments_control(maxit = 0L),
    init = 'zero')
  start = coef(zero); start[] = 0; start[c('sd_re_selection', 'sd_re_outcome')] = 1
  explicit = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', start = start,
    control = arguments_control(maxit = 0L))
  expect_identical(coef(zero), start)
  expect_identical(zero$loglik, explicit$loglik)
})

test_that('panel_control() rejects the fit settings and a carried setting warns', {
  for (name in c('na_action', 'init', 'vcov'))
    expect_error(do.call(panel_control, stats::setNames(list('x'), name)),
      paste(name, 'is now an argument of panel_select()'), fixed = TRUE)
  expect_error(panel_control(vcov = 'opg', init = 'zero'),
    'vcov, init are now arguments of panel_select()', fixed = TRUE)
  data = arguments_data()
  fit = panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = arguments_control(),
    init = 'pooled', vcov = 'observed')
  expect_warning(again <- panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = fit$control),
    paste('control carries init = "pooled", vcov = "observed", but these are panel_select()',
      'arguments; the fit uses init = "default", vcov = "opg"'), fixed = TRUE)
  expect_identical(again$control$vcov, 'opg')
  expect_identical(again$control$init, 'default')
  expect_no_warning(panel_select(y ~ x, d ~ w, data, 'id', 'poisson', control = fit$control,
    init = 'pooled', vcov = 'observed'))
})

test_that('update() keeps the fit settings and accepts them as overrides', {
  data = arguments_data()
  # A missing regressor of a selected row is required data.
  gappy = data; gappy$x[which(data$d == 1L)[2L]] = NA
  fit = panel_select(y ~ x, d ~ w, gappy, 'id', 'poisson', control = arguments_control(),
    init = 'pooled', na_action = 'exclude')
  kept = update(fit)
  expect_identical(kept$control[c('na_action', 'init', 'vcov')],
    list(na_action = 'exclude', init = 'pooled', vcov = 'opg'))
  recipe = update(fit, evaluate = FALSE)
  expect_identical(c(recipe$na_action, recipe$init, recipe$vcov), c('exclude', 'pooled', 'opg'))
  expect_false(any(c('na_action', 'init', 'vcov') %in% names(recipe$control)))
  expect_no_warning(changed <- update(fit, vcov = 'observed', na_action = 'omit'))
  expect_identical(changed$control[c('na_action', 'vcov')],
    list(na_action = 'omit', vcov = 'observed'))
  expect_s3_class(changed$prepared$na.action, 'omit')
  expect_error(update(fit, family = NULL), 'family is required in update()', fixed = TRUE)
  expect_error(update(fit, init = 'nothing', evaluate = FALSE), 'init must be one of')
})

test_that('the default iteration limit is 500 and maxLik loads only when used', {
  expect_identical(panel_control()$maxit, 500L)
  expect_identical(eval(formals(panel_control)$maxit), 500L)
  expect_error(maxlik_namespace(installed = FALSE),
    'optimizer = "maxLik" needs the maxLik package; install it or use optimizer = "optim"',
    fixed = TRUE)
  expect_false(grepl('maxLik', utils::packageDescription('PanelSelect')$Imports))
  expect_match(utils::packageDescription('PanelSelect')$Suggests, 'maxLik')
})

test_that('the default iteration limit lets a Nelder-Mead regression converge', {
  # Ordinary Poisson regression with 13 coefficients needs about 250 Nelder-Mead
  # function evaluations: more than the earlier default of 200.
  local_test_seed(1601)
  n = 300L; k = 12L
  x = matrix(stats::rnorm(n*k), n, k, dimnames = list(NULL, paste0('x', seq_len(k))))
  data = data.frame(id = seq_len(n), x)
  data$count = stats::rpois(n, exp(.2 + drop(x %*% seq(.3, -.3, length.out = k))))
  formula = stats::as.formula(paste('count ~', paste(colnames(x), collapse = ' + ')))
  fit = panel_select(formula, NULL, data, 'id', 'poisson', fixed = c(sd_re_outcome = 0),
    control = panel_control(method = 'Nelder-Mead'))
  expect_identical(fit$optimizer$status, 'converged')
  expect_warning(short <- panel_select(formula, NULL, data, 'id', 'poisson',
    fixed = c(sd_re_outcome = 0), control = panel_control(method = 'Nelder-Mead', maxit = 200L)),
    'did not report convergence')
  expect_identical(short$optimizer$status, 'not_converged')
})

test_that('family binary is the probit outcome', {
  data = arguments_data()
  probit = panel_select(b ~ x, NULL, data, 'id', 'probit', control = arguments_control())
  binary = panel_select(b ~ x, NULL, data, 'id', 'binary', control = arguments_control())
  expect_identical(binary$prepared$family, 'probit')
  expect_identical(coef(binary), coef(probit))
  expect_identical(binary$loglik, probit$loglik)
  # The alias also takes the routes that name the probit family: moving-cut GL
  # for selected binary panels.
  small = data.frame(id = rep(1:3, each = 2), d = c(1, 1, 1, 0, 0, 0), y = c(1, 1, 0, NA, NA, NA))
  theta = c('selection.(Intercept)' = 0, 'outcome.(Intercept)' = 0, sd_re_selection = 1,
    sd_re_outcome = 1, cor_re = 0, cor_error = 0)
  moving_cut = function(family) panel_select(y ~ 1, d ~ 1, small, 'id', family, fixed = theta,
    control = panel_control(quadrature = 'gl'))
  alias = moving_cut('binary')
  expect_identical(alias$integration$method, 'split_legendre')
  expect_identical(alias$loglik_by_unit, moving_cut('probit')$loglik_by_unit)
})

test_that('vcov() accepts vcov = as a synonym of type =', {
  data = arguments_data()
  fit = add_information(panel_select(y ~ x, d ~ w, data, 'id', 'poisson',
    control = arguments_control(), init = 'pooled'))
  expect_identical(vcov(fit, vcov = 'observed'), vcov(fit, type = 'observed'))
  expect_identical(vcov(fit, vcov = 'sandwich'), vcov(fit, 'sandwich'))
  expect_identical(vcov(fit, type = 'opg', vcov = 'opg'), vcov(fit))
  expect_error(vcov(fit, type = 'opg', vcov = 'observed'), 'give one of them')
})
