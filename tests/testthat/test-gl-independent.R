# Independent scalar integrals after collapsing the jointly normal singleton.
# No production contribution, CDF, derivative, or quadrature is used below.
collapsed_singleton_loglik = function(theta, prepared) {
  parameter = function(name) if (name %in% names(theta)) unname(theta[name]) else 0
  b = drop(prepared$x %*% theta[startsWith(names(theta),'outcome.')]) + prepared$offset_outcome
  a = if (prepared$has_selection)
    drop(prepared$w %*% theta[startsWith(names(theta),'selection.')]) + prepared$offset_selection else 0
  lambda = parameter('sd_re_outcome'); delta = parameter('sd_re_selection')
  sigma = parameter('sd_error'); rho = parameter('cor_re'); tau = parameter('cor_error')
  vd = 1+delta^2
  if (prepared$has_selection && prepared$d==0) return(pnorm(-a/sqrt(vd),log.p=TRUE))
  vy = lambda^2 + if (prepared$family=='probit') 1 else sigma^2
  covariance = delta*lambda*rho + if (prepared$family=='probit') tau else sigma*tau
  integrate_scalar = function(f, lower=-Inf, upper=Inf) {
    value = integrate(f,lower,upper,rel.tol=2e-12,abs.tol=0,subdivisions=500L)
    stopifnot(value$message=='OK',value$value>0,value$abs.error/value$value<1e-10)
    value$value
  }
  y = prepared$y
  if (prepared$family=='probit') {
    q = 2*y-1
    if (!prepared$has_selection) return(pnorm(q*b/sqrt(vy),log.p=TRUE))
    correlation = covariance/sqrt(vd*vy)
    return(log(integrate_scalar(function(z) dnorm(z)*
      pnorm((q*b/sqrt(vy)+q*correlation*z)/sqrt(1-correlation^2)),lower=-a/sqrt(vd))))
  }
  if (prepared$family=='gaussian') {
    value = dnorm(y,b,sqrt(vy),log=TRUE)
    return(value+pnorm((a+covariance*(y-b)/vy)/sqrt(vd-covariance^2/vy),log.p=TRUE))
  }
  if(vy==0) return(dpois(y,exp(b),log=TRUE)+
    if(prepared$has_selection) pnorm(a/sqrt(vd),log.p=TRUE) else 0)
  logmass = function(z) {
    index = b+sqrt(vy)*z
    mu = exp(pmin(index,700))
    value = dpois(y,mu,log=TRUE)+dnorm(z,log=TRUE)
    value[index>700] = -Inf
    if (prepared$has_selection) value = value+
      pnorm((a+covariance*z/sqrt(vy))/sqrt(vd-covariance^2/vy),log.p=TRUE)
    value
  }
  log(integrate_scalar(function(z) exp(logmass(z))))
}

independent_singleton_fixture = function(model, sign=1, y=2) {
  cross = model %in% c('probit_probit','probit_linear','probit_pln')
  selected = cross || startsWith(model,'selection_')
  family = if (cross) sub('probit_','',model) else sub('^selection_|_re$','',model)
  if (family=='linear') family = 'gaussian'
  dat = data.frame(id=1L,d=1,y=if(family=='probit') 1 else if(family=='gaussian') .4 else y,
    x=.3,w=-.2,oy=.11,od=-.07)
  prepared = prepare_data(y~x+offset(oy),if(selected) d~w+offset(od) else NULL,
    dat,if(cross) NULL else 'id',family)
  info = parameter_layout(prepared)
  theta = c('selection.(Intercept)'=-.15,selection.w=.3,
    'outcome.(Intercept)'=.25,outcome.x=-.2,sd_re_selection=.45,sd_re_outcome=.6,
    sd_error=.5,cor_re=-.35*sign,cor_error=.25*sign)[info$name]
  list(theta=theta,prepared=prepared)
}

expect_collapsed_singleton = function(model, sign) {
  fixture = independent_singleton_fixture(model,sign)
  control = staged_test_budget(panel_control(quadrature = 'gl',nodes = 15L,
    observation_quadrature = 'gl',observation_nodes = 15L,
    quadrature_options = list(reference_nodes=20L),
    observation_options = list(reference_nodes=20L)))
  if (uses_staged_integration(fixture$prepared,control,numeric()))
    control$integration_mesh = prepare_staged_integration(fixture$theta,fixture$prepared,control,numeric())
  actual = evaluate_likelihood(fixture$theta,fixture$prepared,control,numeric(),score=TRUE)
  reference = collapsed_singleton_loglik(fixture$theta,fixture$prepared)
  expect_lt(abs(actual$loglik-reference),1e-8,label=paste(model,sign,'log error'))
  for (step in c(3e-5,1e-5)) {
    derivative = core_difference(function(theta) collapsed_singleton_loglik(theta,fixture$prepared),
      fixture$theta,step)
    expect_lt(max(abs(drop(actual$score)-derivative)/(1+abs(derivative))),1e-6,
      label=paste(model,sign,step,'scaled score error'))
  }
}

expect_count_singleton = function(model, y) {
  f = independent_singleton_fixture(model,y=y)
  control = staged_test_budget(panel_control(quadrature = 'gl',observation_quadrature = 'gl'))
  control$integration_mesh = prepare_staged_integration(f$theta,f$prepared,control,numeric())
  actual = evaluate_likelihood(f$theta,f$prepared,control,score=FALSE)
  expect_lt(abs(actual$loglik-collapsed_singleton_loglik(f$theta,f$prepared)),1e-8,
    label=paste(model,y,'log error'))
}

test_that('all ten GL engines match independent collapsed singleton values and scores', {
  models = c(core_models,'probit_probit','probit_linear','probit_pln')
  for (model in models) for (sign in c(-1,1)) expect_collapsed_singleton(model,sign)
})

test_that('zero and high count singleton references remain distinct with offsets', {
  for (model in c('poisson_re','pln_re','selection_poisson','selection_pln','probit_pln'))
    for (y in c(0,30)) expect_count_singleton(model,y)
})

# Multi-row count models against the frozen refined references: independent
# base-R likelihoods and scores with 31 outer and 25 inner GH nodes (see
# fixtures/README.md), which GL reproduces here to about 1e-13. Tolerances are
# the staged rules' own default targets, panel_loglik_tol = 1e-7 and panel_score_tol = 1e-5.
# Coupled selected PLN isolates each level: GL outer cells with the reference
# GH25 observation rule, and the reference GH31 outer rule with GL row meshes.
test_that('multi-row GL count models match the frozen refined GH references', {
  control = function(outer, inner) staged_test_budget(panel_control(quadrature = outer,
    nodes = switch(outer, gl = 15L, gh = 31L),
    quadrature_options = if (outer == 'gl') list(reference_nodes = 20L) else list(),
    observation_quadrature = inner, observation_nodes = switch(inner, gl = 15L, gh = 25L),
    observation_options = if (inner == 'gl') list(reference_nodes = 20L) else list()))
  cases = list(c('poisson_re', 'gl', 'gh'), c('pln_re', 'gl', 'gl'),
    c('selection_poisson', 'gl', 'gh'), c('selection_pln', 'gl', 'gh'),
    c('selection_pln', 'gh', 'gl'))
  for (case in cases) {
    fixture = core_fixture(case[1L])
    reference = core_reference(fixture, core_references$refined[[case[1L]]])
    settings = control(case[2L], case[3L])
    settings$integration_mesh = prepare_staged_integration(fixture$theta,
      fixture$prepared, settings, numeric())
    actual = evaluate_model(fixture$theta, fixture$prepared, settings)
    label = paste(case, collapse = ' ')
    expect_lt(abs(sum(actual$loglik)-reference$loglik), 1e-7, label = paste(label, 'log error'))
    expect_lt(max(abs(colSums(actual$score)-reference$score)), 1e-5,
      label = paste(label, 'score error'))
  }
})
