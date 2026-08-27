LL_probitRE_linearRE = function(par,y,d,x,w,group,H=20,rho_off=FALSE, tau_off=FALSE, verbose=1){
  if(tau_off) par = c(par, artanh_tau=0)
  if(rho_off) par = c(par[-length(par)], artanh_rho=0, par[length(par)])
  if(length(par) != ncol(x)+ncol(w)+5){
    print(par)
    stop("Number of parameters incorrect")
  }
  alpha = par[1:ncol(w)]
  beta = par[ncol(w)+1:ncol(x)]
  delta = exp(par[length(par)-4])
  lambda = exp(par[length(par)-3])
  sigma = exp(par[length(par)-2])
  rho = 1 - 2/(exp(par[length(par)-1])+1)
  tau = 1 - 2/(exp(par[length(par)])+1)

  wa = as.vector(w %*% alpha)
  xb = as.vector(x %*% beta)
  rule = gauss.quad(H, "hermite")
  omega = rule$weights
  z = rule$nodes

  Li = rep(0, length(group)-1)
  for(h in 1:H){
    for(k in 1:H){
      A = wa + sqrt(2)*delta*z[k]
      B = xb + sqrt(2)*lambda*(rho*z[k]+sqrt(1-rho^2)*z[h])
      C = (y-B)/sigma
      C[d==0] = 0 # convert NA to zero to avoid NAs
      D = (A+tau*C)/sqrt(1-tau^2)
      # Use log version to avoid numerical issues while taking the product
      log_E = pnorm(D, log.p=TRUE) + dnorm(C, log=TRUE) - log(sigma)
      log_G = d*log_E + (1-d)*pnorm(-A, log.p=TRUE)
      prod = exp(groupSum(log_G, group))
      Li = Li + omega[h]*omega[k]*prod
    }
  }
  Li = pmax(Li/pi, 1e-100) # in case that some Li=0
  LL = sum(log(Li))
  if(verbose>=2){
    writeLines(paste("==== Function call ", panel.select.env$iter, ": LL=",round(LL,digits=5)," =====", sep=""))
    print(round(par,digits=3))
  }
  addIter()
  if(is.na(LL) || !is.finite(LL)){
    if(verbose>=2) writeLines("NA or infinite likelihood, will try others")
    LL = -1e300
  }
  return (LL)
}


Gradient_probitRE_linearRE = function(par,y,d,x,w,group,H=20,rho_off=FALSE, tau_off=FALSE, verbose=1,variance=FALSE){
  if(tau_off) par = c(par, artanh_tau=0)
  if(rho_off) par = c(par[-length(par)], artanh_rho=0, par[length(par)])
  if(length(par) != ncol(x)+ncol(w)+5){
    print(par)
    stop("Number of parameters incorrect")
  }
  alpha = par[1:ncol(w)]
  beta = par[ncol(w)+1:ncol(x)]
  delta = exp(par[length(par)-4])
  lambda = exp(par[length(par)-3])
  sigma = exp(par[length(par)-2])
  rho = 1 - 2/(exp(par[length(par)-1])+1)
  tau = 1 - 2/(exp(par[length(par)])+1)

  wa = as.vector(w %*% alpha)
  xb = as.vector(x %*% beta)
  rule = gauss.quad(H, "hermite")
  omega = rule$weights
  z = rule$nodes

  n = length(group)-1
  Li = rep(0, n)
  dL = matrix(0, n, length(par))
  dG_alpha = matrix(0, n, length(alpha) + 2)
  dG_beta = matrix(0, n, length(beta) + 3)
  w_ext = cbind(w, log_delta=0, artanh_tau=0)
  x_ext = cbind(x, log_lambda=0, log_sigma=0, artanh_rho=0)

  for(h in 1:H){
    for(k in 1:H){
      A = wa + sqrt(2)*delta*z[k]
      B = xb + sqrt(2)*lambda*(rho*z[k]+sqrt(1-rho^2)*z[h])
      C = (y-B)/sigma
      C[d==0] = 0 # convert NA to zero to avoid NAs
      D = (A+tau*C)/sqrt(1-tau^2)
      log_E = pnorm(D, log.p=TRUE) + dnorm(C, log=TRUE) - log(sigma)
      log_G = d*log_E + (1-d)*pnorm(-A, log.p=TRUE)
      prod = exp(groupSum(log_G, group))
      Li = Li + omega[h]*omega[k]*prod

      w_ext[, ncol(w_ext)-1] = sqrt(2)*z[k]
      w_ext[, ncol(w_ext)] = (C+A*tau)/sqrt(1-tau^2)
      x_ext[, ncol(x_ext)-2] = sqrt(2)*(rho*z[k]+sqrt(1-rho^2)*z[h])
      x_ext[, ncol(x_ext)-1] = C
      x_ext[, ncol(x_ext)] = sqrt(2)*lambda*(z[k]-z[h]*rho/sqrt(1-rho^2))

      F_ = dnorm(D)*dnorm(C)/(sigma*sqrt(1-tau^2))
      E = exp(log_E)
      G = exp(log_G)
      dG_alpha = matVecProd(w_ext, d*F_) + cbind(matVecProd(w_ext[,-ncol(w_ext)], -(1-d)*dnorm(-A)), 0)
      dG_beta = matVecProd(x_ext, -d*(F_*tau-E*C)/sigma)
      dG_beta[, ncol(x_ext)-1] = dG_beta[, ncol(x_ext)-1] - d*E/sigma

      sumT = matVecProdSum(cbind(dG_alpha, dG_beta), numeric(0), 1/G, group)
      dL = dL + matVecProd(sumT, omega[h]*omega[k]*prod)
    }
  }
  Li = pmax(Li/pi, 1e-100) # in case that some Li=0
  dL = matVecProd(dL, 1/(pi*Li))
  # rearrange columns
  # colnames(dL) = c(colnames(w_ext), colnames(x_ext))
  dL = dL[, c(1:ncol(w), ncol(w_ext)+1:ncol(x), ncol(w)+1, ncol(w_ext)+ncol(x)+1:3, ncol(w)+2)]
  # print(colnames(dL))
  # print(names(par))
  colnames(dL) = names(par)

  # accounting for transformation of parameters
  dL[, 'log_delta'] = delta * dL[, 'log_delta']
  dL[, 'log_lambda'] = lambda * dL[, 'log_lambda']
  dL[, 'log_sigma'] = sigma * dL[, 'log_sigma']
  dL[, 'artanh_rho'] = (2 * exp(rho) / (exp(rho)+1)^2) * dL[, 'artanh_rho']
  dL[, 'artanh_tau'] = (2 * exp(tau) / (exp(tau)+1)^2) * dL[, 'artanh_tau']

  # Remove rho and tau if they are off
  if(rho_off) dL = dL[, -(ncol(dL)-1)]
  if(tau_off) dL = dL[, -ncol(dL)]

  gradient = colSums(dL)

  # num_g = numericGradient(LL_probitRE_linearRE, par, y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off)
  # cat('-------Gradient difference------\n')
  # print(num_g - gradient)

  if(verbose>=3){
    cat("----Gradient:\n")
    print(gradient,digits=3)
  }
  if(any(is.na(gradient) | !is.finite(gradient))) gradient = rep(NA, length(gradient))
  if(variance){
    var = tryCatch( solve(crossprod(dL)), error = function(e){
      cat('BHHH cross-product not invertible: ', e$message, '\n')
      diag(length(par)) * NA
    } )
    return (list(g=gradient, var=var, I=crossprod(dL)))
  }
  return(gradient)
}

#' Prediction Function for probitRE_linearRE Models
#'
#' Generates predictions from a fitted \code{probitRE_linearRE} model.
#'
#' @param model A fitted \code{probitRE_linearRE} model object.
#' @param newdata A data frame containing the covariates used for prediction.
#'
#' @return A list of predicted values.
#'
#' @export
#' @md
predict_probitRE_linearRE = function(model, newdata){
  par = model$estimates[, 1]
  var = model$var
  form_probit = model$form_probit
  form_linear = model$form_linear

  # parse linear formula
  mf = model.frame(form_linear, data=newdata, na.action=NULL, drop.unused.levels=TRUE)
  y = model.response(mf, "numeric")
  x = model.matrix(attr(mf, "terms"), data=mf)

  # parse probit formula
  mf2 = model.frame(form_probit, data=newdata, na.action=NULL, drop.unused.levels=TRUE)
  d = model.response(mf2, "numeric")
  w = model.matrix(attr(mf2, "terms"), data=mf2)

  # separate parameters
  alpha = par[1:ncol(w)]
  beta = par[ncol(w)+1:ncol(x)]
  delta = par['delta']
  lambda = par['lambda']
  sigma = par['sigma']
  rho = ifelse("rho" %in% names(par), par["rho"], 0)
  tau = ifelse("tau" %in% names(par), par["tau"], 0)

  wa = as.vector(w %*% alpha)
  xb = as.vector(x %*% beta)

  respond_prob = pnorm(wa/sqrt(1+delta^2))
  outcome = xb
  # list(respond_prob = respond_prob, outcome = outcome, d=d)

  # SE of outcome
  ix2 = ncol(w)+1:ncol(x)
  gr_outcome = x
  se_outcome = apply(gr_outcome, 1, function(g) drop(sqrt(t(g) %*% var[ix2,ix2] %*% g)))

  # SE of respond_prob
  ix1 = c(1:ncol(w), ncol(w)+ncol(x)+1)
  gr_respond = matVecProd(cbind(w, -delta*wa/(1+delta^2)), respond_prob/sqrt(1+delta^2))
  se_respond = apply(gr_respond, 1, function(g) drop(sqrt(t(g) %*% var[ix1,ix1] %*% g)))

  # Population mean
  pop_mean = mean(outcome)
  gr_mean = colMeans(gr_outcome)
  pop_mean_se = drop(sqrt(t(gr_mean) %*% var[ix2,ix2] %*% gr_mean))

  list(respond_prob = compile(respond_prob, se_respond),
       outcome = compile(outcome, se_outcome),
       pop_mean = compile(pop_mean, pop_mean_se),
       gr_respond = gr_respond, gr_outcome = gr_outcome, d=d, y=y)
}


#' Panel Sample Selection Model for Continuous Outcome
#' @description A panel sample selection model for continuous outcome, with selection at both the individual and individual-time levels. The outcome is observed in the second stage only if the first stage outcome is one.\cr\cr
#' Let \eqn{\boldsymbol{w}_{it}} and \eqn{\boldsymbol{x}_{it}} represent the *row* vectors of covariates in the selection and outcome equations, respectively, with \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}} denoting the corresponding *column* vectors of parameters.\cr\cr
#' First stage (probitRE):
#' \deqn{d_{it}=1(\mathbf{w}_{it} \boldsymbol{\alpha}+\delta u_i+\varepsilon_{it}>0)}{d_it = 1(w_it * \alpha + \delta * u_i +\varepsilon_it > 0)}
#' Second stage (linearRE):
#' \deqn{y_{it} = \mathbf{x}_{it} \boldsymbol{\beta} + \lambda v_i +\sigma \epsilon_{it}}{y_it = x_it * \beta + \gamma * m_i + \lambda * v_i + \sigma * \epsilon_it}
#' Correlation structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}.
#' \eqn{\varepsilon_{it}} and \eqn{\epsilon_{it}} are bivariate normally distributed with a correlation of \eqn{\tau}. \cr\cr
#' w and x can be the same set of variables. Identification can be weak if w are not good predictors of d.
#' @param form_probit Formula for the panel probit model with random effects at the individual level
#' @param form_linear Formula for the panel linear model with random effects at the individual level
#' @param id.name the name of the id column in data
#' @param data Input data, must be a data.frame object
#' @param par Starting values for estimates
#' @param init Initialization method
#' @param method Optimization algorithm. Default is BFGS
#' @param H Number of quadrature points
#' @param rho_off A Boolean value indicating whether to turn off the correlation between the random effects of the probit and linear models. Default is FALSE.
#' @param tau_off A Boolean value indicating whether to turn off the correlation between the error terms of the probit and linear models. Default is FALSE.
#' @param use.optim A Boolean value indicating whether to use optim instead of maxLik. Default is FALSE.
#' @param rho.init Initial value for the correlation between the random effects of the probit and linear models. Default is 0.
#' @param tau.init Initial value for the correlation between the error terms of the probit and linear models. Default is 0.
#' @param verbose A integer indicating how much output to display during the estimation process.
#' * <0 - No ouput
#' * 0 - Basic output (model estimates)
#' * 1 - Limited output, providing likelihood of iterations
#' * 2 - Moderate output, basic ouput + parameter and likelihood on each call
#' * 3 - Extensive output, moderate output + gradient values on each call
#' @return A list containing the results of the estimated model, some of which are inherited from the return of maxLik
#' * estimates: Model estimates with 95% confidence intervals
#' * estimate or par: Point estimates
#' * predict. A list containing the predicted probabilities of responding (respond_prob) and the predicted counterfactual outcome values (outcome), their gradients (gr_respond and gr_outcome), and estimated counterfactual population mean (pop_mean).
#' * variance_type: covariance matrix used to calculate standard errors. Either BHHH or Hessian.
#' * var: covariance matrix
#' * se: standard errors
#' * var_bhhh: BHHH covariance matrix, inverse of the outer product of gradient at the maximum
#' * se_bhhh: BHHH standard errors
#' * gradient: Gradient function at maximum
#' * hessian: Hessian matrix at maximum
#' * gtHg: \eqn{g'H^-1g}, where H^-1 is simply the covariance matrix. A value close to zero (e.g., <1e-3 or 1e-6) indicates good convergence.
#' * LL or maximum: Likelihood
#' * AIC: AIC
#' * BIC: BIC
#' * n_obs: Number of observations
#' * n_par: Number of parameters
#' * time: Time takes to estimate the model
#' * iterations: number of iterations taken to converge
#' * message: Message regarding convergence status.
#'
#' Note that the list inherits all the components in the output of maxLik. See the documentation of maxLik for more details.
#' @md
#' @examples
#' library(PanelSelect)
#' library(MASS)
#' N = 200
#' period = 3
#' obs = N*period
#' rho = 0.5
#' tau = 0
#' set.seed(1)
#'
#' re = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
#' u = rep(re[,1], each=period)
#' v = rep(re[,2], each=period)
#' e = mvrnorm(obs, mu=c(0,0), Sigma=matrix(c(1,tau,tau,1), nrow=2))
#' e1 = e[,1]
#' e2 = e[,2]
#'
#' t = rep(1:period, N)
#' id = rep(1:N, each=period)
#' w = rnorm(obs)
#' z = rnorm(obs)
#' x = rnorm(obs)
#' d = as.numeric(x + w + u + e1 > 0)
#' y = x + w + v + e2
#' y[d==0] = NA
#' dt = data.frame(id, t, y, x, w, z, d)
#'
#' # As N increases, the parameter estimates will be more accurate
#' # Do not turn off tau if you believe it is nonzero
#' m = probitRE_linearRE(d~x+w, y~x+w, 'id', dt, H=10, tau_off=TRUE, verbose=-1)
#' print(m$estimates, digits=4)
#' @export
#' @family PanelSelect
#' @references
#' Bailey, M., & Peng, J. (2025). A Random Effects Model of Non-Ignorable Nonresponse in Panel Survey Data. Available at SSRN <https://www.ssrn.com/abstract=5475626>
probitRE_linearRE = function(form_probit, form_linear, id.name, data=NULL, par=NULL, method='BFGS', rho_off=FALSE, tau_off=FALSE, H=10, init=c('zero', 'unif', 'norm', 'default')[4], rho.init=0, tau.init=0, use.optim=FALSE, verbose=0){
  # 1.1 Sort data by id
  data = data.frame(data)
  data_original = copy(data)
  ord = order(data[, id.name])
  data = data[ord, ]
  group = c(0,cumsum(table(as.integer(factor(data[, id.name])))))

  # 1.1 parse linear formula
  mf = model.frame(form_linear, data=data, na.action=NULL, drop.unused.levels=TRUE)
  y = model.response(mf, "numeric")
  x = model.matrix(attr(mf, "terms"), data=mf)

  # 1.2 parse probit formula
  mf2 = model.frame(form_probit, data=data, na.action=NULL, drop.unused.levels=TRUE)
  d = model.response(mf2, "numeric")
  w = model.matrix(attr(mf2, "terms"), data=mf2)

  # 1.3 Initialize parameters
  est_linear = lm(form_linear, data=data)
  par_linear = coef(summary(est_linear))[,1]
  est_probit = glm(form_probit, data=data, family=binomial(link="probit"))
  par_probit = coef(summary(est_probit))[,1]
  names(par_linear) = paste0('linear.', names(par_linear))
  names(par_probit) = paste0('probit.', names(par_probit))
  par_linear[is.na(par_linear)] = 0
  par_probit[is.na(par_probit)] = 0
  # convert bounded parameters to unboudned, rho to 2*atanh(rho)=ln((1+rho)/(1-rho))
  par = c(par_probit, par_linear, log_delta=0, log_lambda=0, log_sigma=0, artanh_rho=2*atanh(rho.init), artanh_tau=2*atanh(tau.init))
  if(init=='unif') par = par - par + runif(length(par))
  if(init=='norm') par = par - par + rnorm(length(par))
  if(init=='zero') par = par - par
  if(rho_off) par = par[-(length(par)-1)]
  if(tau_off) par = par[-length(par)]
  # print(par)

  # 2. Estimation
  panel.select.env$LL = -.Machine$double.xmax
  panel.select.env$iter = 1
  begin = Sys.time()

  # use maxLik (identical estimate with optim, but more reliable SE)
  if(use.optim){
    res = optim(par=par, fn=LL_probitRE_linearRE, gr=Gradient_probitRE_linearRE, method=method, control=list(fnscale=-1, trace=as.numeric(verbose>=1)), y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off, verbose=verbose, hessian = TRUE)
    res$iterations = res$counts['function']
    res$LL = res$value
  }else{
    res = maxLik(LL_probitRE_linearRE, grad=Gradient_probitRE_linearRE, start=par, y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off, method=method, verbose=verbose, printLevel=verbose+1)
    res$par = res$estimate
    res$LL = res$maximum
  }
  res$n_obs = length(d)

  # 3. Compile results
  # res = getVarSE(res, verbose=verbose)
  gvar = Gradient_probitRE_linearRE(res$par,y,d,x,w,group,H,rho_off=rho_off,tau_off=tau_off,verbose=verbose-1,variance=TRUE)
  if(use.optim) res$gradient = gvar$g # optim has no gradient
  res = getVarSE(res, gvar=gvar, verbose=verbose)

  # res$num_g = numericGradient(LL_probitRE_linearRE,res$par,y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off)
  # cat('-------Gradient difference------\n')
  # print(res$num_g - gvar$g)

  # convert unbounded estimates to bounded estimates
  trans_vars=c(delta='log_delta', lambda='log_lambda', sigma='log_sigma', rho='artanh_rho', tau='artanh_tau')
  trans_types=c('exp', 'exp', 'exp', 'correlation', 'correlation')
  if(rho_off){
    trans_vars = trans_vars[-(length(trans_vars)-1)]
    trans_types = trans_types[-(length(trans_types)-1)]
  }
  if(tau_off){
    trans_vars = trans_vars[-length(trans_vars)]
    trans_types = trans_types[-length(trans_types)]
  }
  res = transCompile(res, trans_vars, trans_types)
  res$form_probit = form_probit
  res$form_linear = form_linear
  res$predict = predict_probitRE_linearRE(res, data_original)

  # Need to estimate probitRE and linearRE models to make the test meaningful
  # res$LR_stat = 2 * ( res$LL - logLik(est_linear) - logLik(est_probit) )
  # res$LR_p = 1 - pchisq(res$LR_stat, 1)
  res$ord = ord
  res$iter = panel.select.env$iter

  if(verbose>=0){
    cat(sprintf('==== Converged after %d calls of likelihood function, LL=%.2f, gtHg=%.6f ****\n', res$iterations, res$LL, res$gtHg))
    # cat(sprintf('LR test of rho=0, chi2(1)=%.3f, p-value=%.4f\n', res$LR_stat, res$LR_p))
    print(res$time <- Sys.time() - begin)
  }
  return (res)
}

