# History predictions and posterior random effects against direct tensor
# integration in reference/math.R and reference/prediction.R (no production
# code). Fits fix every parameter, so the comparisons are at known values.

history_reference_case = function(model) {
  spec = ref_model(model)
  data = data.frame(id = rep(1:3, each = 3),
    x = c(-.6, .2, .9, .4, -.3, .1, -1, .5, .7),
    w = c(.3, -.5, .8, .1, .6, -.2, -.4, .9, .2),
    d = c(1, 0, 1, 1, 1, 0, 0, 1, 1))
  data$y = switch(spec$family, probit = c(1, NA, 0, 1, 0, NA, NA, 1, 1), c(2, NA, 0, 1, 3, NA, NA, 0, 4))
  if (!spec$selection) {
    data$d = NULL
    data$y[is.na(data$y)] = 1
  }
  theta = c('outcome.(Intercept)' = .2, outcome.x = .4, sd_re_outcome = .6)
  if (spec$selection) theta = c('selection.(Intercept)' = .1, selection.w = .5, theta,
    sd_re_selection = .7, cor_re = .4)
  if (spec$family == 'pln') theta = c(theta, sd_error = .5)
  if (spec$selection && spec$family != 'poisson') theta = c(theta, cor_error = -.3)
  par = theta
  names(par)[startsWith(names(par), 'outcome.')] = c('beta1', 'beta2')
  names(par)[startsWith(names(par), 'selection.')] = c('alpha1', 'alpha2')
  panels = lapply(split(data, data$id), function(rows) list(x = cbind(1, rows$x),
    w = cbind(1, rows$w), y = rows$y, d = rows$d))
  list(spec = spec, data = data, theta = theta, par = par, panels = panels)
}

test_that('history predictions equal directly integrated conditional moments', {
  nodes = 24L
  for (model in c('selection_probit', 'selection_pln', 'poisson_re')) {
    case = history_reference_case(model)
    spec = case$spec
    fit = panel_select(y ~ x, if (spec$selection) d ~ w, case$data, 'id', spec$family,
      fixed = case$theta, control = panel_control(quadrature = 'gh', nodes = nodes,
        observation_nodes = 30L))
    targets = data.frame(id = 1:3, x = c(.3, -.2, .5), w = c(-.2, .4, .1))
    # Reference keys and the prediction targets they check.
    types = c(response = 'response', if (spec$selection) c(selection = 'selection_probability',
      joint = 'joint', selected = 'selected_mean'))
    for (k in 1:3) {
      expected = ref_prediction_integral(model, case$par, c(1, targets$x[k]), c(1, targets$w[k]),
        case$panels[[k]], outer = nodes, inner = 30L)
      for (key in names(types)) {
        actual = predict(fit, targets[k, ], type = types[[key]], conditioning = 'history',
          history = case$data)
        expect_equal(unname(actual), unname(expected[[key]]), tolerance = 1e-8,
          info = paste(model, k, key))
      }
    }
  }
})

test_that('random effects equal directly integrated posterior moments', {
  nodes = 24L
  for (model in c('selection_probit', 'selection_pln', 'poisson_re')) {
    case = history_reference_case(model)
    spec = case$spec
    fit = panel_select(y ~ x, if (spec$selection) d ~ w, case$data, 'id', spec$family,
      fixed = case$theta, control = panel_control(quadrature = 'gh', nodes = nodes,
        observation_nodes = 30L))
    effects = random_effects(fit, case$data)
    expect_identical(effects$id, 1:3)
    for (k in 1:3) {
      expected = ref_posterior_effects(model, case$par, case$panels[[k]], outer = nodes,
        inner = 30L)
      # Columns after the ID follow selection mean/SD (selection models), then outcome mean/SD.
      expect_equal(unname(unlist(effects[k, -1L])), unname(expected), tolerance = 1e-8,
        info = paste(model, k))
    }
  }
})
