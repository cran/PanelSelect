# Multiwave normal-probability references retain all selection events and only
# observed outcomes. Miwa refinement is corroboration, not certification.
multiwave_mesh = function(fixture,control,fixed=numeric()) {
  control$integration_mesh=prepare_staged_integration(fixture$theta,fixture$prepared,control,fixed)
  control
}
multiwave_normal_reference = function(theta, prepared, steps=1024L) {
  parameter = function(name) if(name %in% names(theta)) unname(theta[name]) else 0
  a = drop(prepared$w %*% theta[startsWith(names(theta),'selection.')])+prepared$offset_selection
  b = drop(prepared$x %*% theta[startsWith(names(theta),'outcome.')])+prepared$offset_outcome
  n = length(a); yes = which(prepared$d==1); r = length(yes); qd = 2*prepared$d-1
  delta=parameter('sd_re_selection'); lambda=parameter('sd_re_outcome')
  sigma=parameter('sd_error'); rho=parameter('cor_re'); tau=parameter('cor_error')
  vd = matrix(delta^2,n,n)+diag(n)
  cross = matrix(delta*lambda*rho,n,r)+
    if(prepared$family=='gaussian') sigma*tau*diag(n)[,yes,drop=FALSE] else tau*diag(n)[,yes,drop=FALSE]
  if(prepared$family=='probit') {
    vy = matrix(lambda^2,r,r)+diag(r)
    signs = c(qd,2*prepared$y[yes]-1)
    mean = c(a,b[yes])*signs
    covariance = rbind(cbind(vd,cross),cbind(t(cross),vy))*tcrossprod(signs)
    return(log(as.numeric(mvtnorm::pmvnorm(lower=-mean,upper=Inf,
      sigma=covariance,algorithm=mvtnorm::Miwa(steps=steps)))))
  }
  vy = matrix(lambda^2,r,r)+diag(sigma^2,r)
  residual = prepared$y[yes]-b[yes]
  solved = solve(vy,residual)
  marginal = -(r*log(2*pi)+as.numeric(determinant(vy,logarithm=TRUE)$modulus)+sum(residual*solved))/2
  mean = (a+drop(cross %*% solved))*qd
  covariance = (vd-cross %*% solve(vy,t(cross)))*tcrossprod(qd)
  marginal+log(as.numeric(mvtnorm::pmvnorm(lower=-mean,upper=Inf,
    sigma=covariance,algorithm=mvtnorm::Miwa(steps=steps))))
}

test_that('multiwave coupled binary and Gaussian values match external normal probabilities', {
  skip_if_not_installed('mvtnorm')
  for(model in c('selection_probit','selection_gaussian')) for(sign in c(-1,1)) {
    f=core_fixture(model,sign)
    f$prepared$offset_outcome=c(.07,-.08,.13); f$prepared$offset_selection=c(-.11,.02,.06)
    references=vapply(c(256L,512L,1024L,2048L,4096L),function(n)
      multiwave_normal_reference(f$theta,f$prepared,n),numeric(1))
    # Coarse Miwa resolutions can differ by 5e-8 here; qualification uses the
    # final two resolutions and retains the same 1e-9 reference-spread budget.
    expect_lt(diff(range(tail(references,2))),1e-9)
    control=multiwave_mesh(f,staged_test_budget(panel_control(quadrature = 'gl',
      nodes = 15L, quadrature_options = list(reference_nodes=20L))))
    value=evaluate_likelihood(f$theta,f$prepared,control)
    expect_lt(abs(value$loglik-tail(references,1)),1e-8,label=paste(model,sign,'log error'))
  }
})

combinations=expand.grid(outer=c('gh','agh','gl'),inner='gl',
  stringsAsFactors=FALSE)
for(k in seq_len(nrow(combinations))) {
  test_that(paste('selected PLN finite-stage scores use',combinations$outer[k],
    'outer and',combinations$inner[k],'inner integration'), {
    f=core_fixture('selection_pln')
    info=parameter_layout(f$prepared)
    f$prepared$offset_outcome=c(.07,-.08,.13); f$prepared$offset_selection=c(-.11,.02,.06)
    fixed=if(k %% 2L) c(cor_error=0) else numeric()
    if(length(fixed)) f$theta[names(fixed)]=fixed
    if(k %in% c(3L,7L)) f$theta['cor_re']=0
    outer=combinations$outer[k]; inner=combinations$inner[k]
    # These finite-objective derivative checks retain all numerical settings,
    # with generous time budgets so allocation speed does not define accuracy.
    control=staged_test_budget(panel_control(quadrature = outer,
      nodes = if(outer=='gl') 15L else 7L,
      observation_quadrature = inner,observation_nodes = 15L,
      quadrature_options = if(outer=='gl') list(reference_nodes=20L) else list(),
      observation_options = list(reference_nodes=20L)))
    control=multiwave_mesh(f,control,fixed)
    saved=serialize(control$integration_mesh,NULL)
    actual=evaluate_likelihood(f$theta,f$prepared,control,fixed,TRUE)
    frozen=control; frozen$.staged_preparing=TRUE
    info$fixed=info$name %in% names(fixed)
    for(step in c(3e-5,1e-5)) {
      free=!names(f$theta) %in% names(fixed)
      # A fixed coordinate is part of the model restriction. Perturb only free
      # coordinates, as the optimizer does, keeping each retained mesh fixed.
      numeric=core_difference(function(value) {
        theta=f$theta;theta[free]=value
        sum(evaluate_likelihood(theta,f$prepared,frozen,fixed,FALSE)$loglik)
      },f$theta[free],step)
      expect_lt(max(abs(colSums(actual$score)[free]-numeric)/(1+abs(numeric))),2e-6,
        label=paste(outer,inner,step,'natural score error'))
      eta=parameter_to_optimizer(f$theta,info)
      numeric_eta=core_difference(function(x) sum(evaluate_likelihood(
        parameter_from_optimizer(x,f$theta,info),f$prepared,frozen,fixed,FALSE)$loglik),eta,step)
      gradient=colSums(actual$score)[free]*parameter_jacobian(f$theta,info)[free]
      expect_lt(max(abs(gradient-numeric_eta)/(1+abs(numeric_eta))),2e-6,
        label=paste(outer,inner,step,'optimizer score error'))
    }
    expect_identical(saved,serialize(control$integration_mesh,NULL))
    again=evaluate_likelihood(f$theta,f$prepared,control,fixed,FALSE)
    expect_equal(again$loglik,actual$loglik,tolerance=1e-12)
  })
}
