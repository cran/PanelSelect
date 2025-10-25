LL_probitRE_probitRE = function(par,y,d,x,w,group,H=20,rho_off=FALSE, tau_off=FALSE, verbose=1){
  if(tau_off) par = c(par, artanh_tau=0)
  if(rho_off) par = c(par[-length(par)], artanh_rho=0, par[length(par)])
  if(length(par) != ncol(x)+ncol(w)+4){
    print(par)
    stop("Number of parameters incorrect")
  }
  alpha = par[1:ncol(w)]
  beta = par[ncol(w)+1:ncol(x)]
  delta = exp(par[length(par)-3])
  lambda = exp(par[length(par)-2])
  rho = 1 - 2/(exp(par[length(par)-1])+1)
  tau = 1 - 2/(exp(par[length(par)])+1)

  wa = as.vector(w %*% alpha)
  xb = as.vector(x %*% beta)
  sign_y = 2*y[d==1] - 1 # sign of 2y-1 when d=1
  rule = gauss.quad(H, "hermite")
  omega = rule$weights
  z = rule$nodes

  Li = rep(0, length(group)-1)
  for(h in 1:H){
    for(k in 1:H){
      A = wa + sqrt(2)*delta*z[k]
      B = xb + sqrt(2)*lambda*(rho*z[k]+sqrt(1-rho^2)*z[h])
      log_P = rep(0, length(d))
      log_P[d==0] = pnorm(-A[d==0], log.p=TRUE)
      log_P[d==1] = log(pmax(pbivnorm(A[d==1], sign_y*B[d==1], sign_y*tau), 1e-100))
      # if(sum(is.na(log_P))){
      #   print(par, digits=3)
      #   # pbvnorm is faster but often produces NA when value is close to 0or 1
      #   Add "@importFrom pbv pbvnorm" to PanelSelect.R if enabled
      #   p1 = pbvnorm(A[d==1], sign_y*B[d==1], sign_y*tau)
      #   p2 = pbivnorm(A[d==1], sign_y*B[d==1], sign_y*tau)
      #   print(cbind(p1=as.vector(p1), p2=as.vector(p2), A=A[d==1], B=B[d==1], tau=tau)[is.na(p1), ], digits=3)
      #   # stop('NA in log_P')
      # }
      prod = exp(groupSum(log_P, group))
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


Gradient_probitRE_probitRE = function(par,y,d,x,w,group,H=20,rho_off=FALSE, tau_off=FALSE, verbose=1,variance=FALSE){
  if(tau_off) par = c(par, artanh_tau=0)
  if(rho_off) par = c(par[-length(par)], artanh_rho=0, par[length(par)])
  if(length(par) != ncol(x)+ncol(w)+4){
    print(par)
    stop("Number of parameters incorrect")
  }
  alpha = par[1:ncol(w)]
  beta = par[ncol(w)+1:ncol(x)]
  delta = exp(par[length(par)-3])
  lambda = exp(par[length(par)-2])
  rho = 1 - 2/(exp(par[length(par)-1])+1)
  tau = 1 - 2/(exp(par[length(par)])+1)

  wa = as.vector(w %*% alpha)
  xb = as.vector(x %*% beta)
  sign_y = 2*y[d==1] - 1 # sign of 2y-1 when d=1
  rule = gauss.quad(H, "hermite")
  omega = rule$weights
  z = rule$nodes

  n = length(group)-1
  Li = rep(0, n)
  dL = matrix(0, n, length(par))
  dG_alpha = matrix(0, n, length(alpha) + 1)
  dG_beta = matrix(0, n, length(beta) + 2)
  w_ext = cbind(w, log_delta=0)
  x_ext = cbind(x, log_lambda=0, artanh_rho=0)

  for(h in 1:H){
    for(k in 1:H){
      A = wa + sqrt(2)*delta*z[k]
      B = xb + sqrt(2)*lambda*(rho*z[k]+sqrt(1-rho^2)*z[h])
      A0 = A[d==0]
      A1 = A[d==1]
      B0 = B[d==0]
      B1 = B[d==1]

      log_P = rep(0, length(d))
      log_P[d==0] = pnorm(-A0, log.p=TRUE)
      log_P[d==1] = log(pmax(pbivnorm(A[d==1], sign_y*B[d==1], sign_y*tau), 1e-100))
      prod = exp(groupSum(log_P, group))
      Li = Li + omega[h]*omega[k]*prod

      w_ext[, ncol(w_ext)] = sqrt(2)*z[k]
      x_ext[, ncol(x_ext)-1] = sqrt(2)*(rho*z[k]+sqrt(1-rho^2)*z[h])
      x_ext[, ncol(x_ext)] = sqrt(2)*lambda*(z[k]-z[h]*rho/sqrt(1-rho^2))

      C = D = E = rep(0, length(d))
      C[d==0] = - dnorm(A0)
      C[d==1] = dnorm(A1)*pnorm(sign_y*(B1-tau*A1)/sqrt(1-tau^2))
      D[d==1] = sign_y*dnorm(sign_y*B1)*pnorm((A1-tau*B1)/sqrt(1-tau^2))
      E[d==1] = sign_y*dbvnorm(A1, sign_y*B1, sign_y*tau)

      dG_alpha = matVecProd(w_ext, C)
      dG_beta = matVecProd(x_ext, D)
      sumT = matVecProdSum(cbind(dG_alpha, dG_beta, E), numeric(0), 1/exp(log_P), group)
      dL = dL + matVecProd(sumT, omega[h]*omega[k]*prod)
    }
  }
  Li = pmax(Li/pi, 1e-100) # in case that some Li=0
  dL = matVecProd(dL, 1/(pi*Li))
  # rearrange columns
  dL = dL[, c(1:ncol(w), ncol(w_ext)+1:ncol(x), ncol(w)+1, ncol(w_ext)+ncol(x)+1:2, ncol(dL))]
  colnames(dL) = names(par)

  # accounting for transformation of parameters
  dL[, 'log_delta'] = delta * dL[, 'log_delta']
  dL[, 'log_lambda'] = lambda * dL[, 'log_lambda']
  dL[, 'artanh_rho'] = (2 * exp(rho) / (exp(rho)+1)^2) * dL[, 'artanh_rho']
  dL[, 'artanh_tau'] = (2 * exp(tau) / (exp(tau)+1)^2) * dL[, 'artanh_tau']

  # Remove rho and tau if they are off
  if(rho_off) dL = dL[, -(ncol(dL)-1)]
  if(tau_off) dL = dL[, -ncol(dL)]

  gradient = colSums(dL)

  # num_g = numericGradient(LL_probitRE_probitRE, par, y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off)
  # cat('-------gradient, num_g, ratio------\n')
  # print(cbind(gradient, g=as.vector(num_g), ratio=gradient/as.vector(num_g)), digits=3)

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

predict_probitRE_probitRE = function(par,var,probit1,probit2,data){
  # parse outcome formula
  mf = model.frame(probit2, data=data, na.action=NULL, drop.unused.levels=TRUE)
  y = model.response(mf, "numeric")
  x = model.matrix(attr(mf, "terms"), data=mf)

  # parse selection formula
  mf2 = model.frame(probit1, data=data, na.action=NULL, drop.unused.levels=TRUE)
  d = model.response(mf2, "numeric")
  w = model.matrix(attr(mf2, "terms"), data=mf2)

  # separate parameters
  alpha = par[1:ncol(w)]
  beta = par[ncol(w)+1:ncol(x)]
  delta = par['delta']
  lambda = par['lambda']
  rho = ifelse("rho" %in% names(par), par["rho"], 0)
  tau = ifelse("tau" %in% names(par), par["tau"], 0)

  wa = as.vector(w %*% alpha)
  xb = as.vector(x %*% beta)

  respond_prob = pnorm(wa/sqrt(1+delta^2))
  outcome_prob = pnorm(xb/sqrt(1+lambda^2))

  # SE of outcome_prob
  ix2 = c(ncol(w)+1:ncol(x), ncol(w)+ncol(x)+2)
  gr_outcome = matVecProd(cbind(x, -lambda*xb/(1+lambda^2)), outcome_prob/sqrt(1+lambda^2))
  se_outcome = apply(gr_outcome, 1, function(g) drop(sqrt(t(g) %*% var[ix2,ix2] %*% g)))

  # SE of respond_prob
  ix1 = c(1:ncol(w), ncol(w)+ncol(x)+1)
  gr_respond = matVecProd(cbind(w, -delta*wa/(1+delta^2)), respond_prob/sqrt(1+delta^2))
  se_respond = apply(gr_respond, 1, function(g) drop(sqrt(t(g) %*% var[ix1,ix1] %*% g)))

  # Population mean
  pop_mean = mean(outcome_prob)
  gr_mean = colMeans(gr_outcome)
  pop_mean_se = drop(sqrt(t(gr_mean) %*% var[ix2,ix2] %*% gr_mean))

  list(respond_prob = compile(respond_prob, se_respond),
       outcome_prob = compile(outcome_prob, se_outcome),
       # CI and p-value may be unreliable as pop_mean lies in [0,1] and unlikely normal
       pop_mean = compile(pop_mean, pop_mean_se),
       gr_respond = gr_respond, gr_outcome = gr_outcome, d=d, y=y)
}

#' Panel Sample Selection Model for Binary Outcome
#' @description A panel sample selection model for binary outcome, with selection at both the individual and individual-time levels. The outcome is observed in the second stage only if the first stage outcome is one.\cr\cr
#' Let \eqn{\mathbf{w}_{it}} and \eqn{\mathbf{x}_{it}} represent the *row* vectors of covariates in the selection and outcome equations, respectively, with \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}} denoting the corresponding *column* vectors of parameters.\cr\cr
#' First stage (probitRE):
#' \deqn{d_{it}=1(\mathbf{w}_{it} \boldsymbol{\alpha}+\delta u_i+\varepsilon_{it}>0)}{d_it = 1(w_it * \alpha + \delta * u_i +\varepsilon_it > 0)}
#' Second stage (probitRE):
#' \deqn{y_{it} = 1(\mathbf{x}_{it} \boldsymbol{\beta} + \lambda v_i +\epsilon_{it}>0)}{y_it = 1(x_it * \beta + \gamma * m_i + \lambda * v_i + \epsilon_it > 0)}
#' Correlation structure:
#' \eqn{u_i} and \eqn{v_i} are bivariate normally distributed with a correlation of \eqn{\rho}.
#' \eqn{\varepsilon_{it}} and \eqn{\epsilon_{it}} are bivariate normally distributed with a correlation of \eqn{\tau}. \cr\cr
#' w and x can be the same set of variables. Identification can be weak if w are not good predictors of d.
#' @param probit1 Formula for the first-stage probit model with random effects at the individual level
#' @param probit2 Formula for the second-stage probit model with random effects at the individual level
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
#' * predict. A list containing the predicted probabilities of responding (respond_prob) and the predicted counterfactual outcome values (outcome_prob), their gradients (gr_respond and gr_outcome), and estimated counterfactual population mean (pop_mean).
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
#' N = 150
#' period = 5
#' obs = N*period
#' rho = 0.5
#' tau = 0
#' set.seed(100)
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
#' y = as.numeric(x + w + v + e2 > 0)
#' y[d==0] = NA
#' dt = data.frame(id, t, y, x, w, z, d)
#'
#' # As N increases, the parameter estimates will be more accurate
#' m = probitRE_probitRE(d~x+w, y~x+w, 'id', dt, H=10, verbose=-1)
#' print(m$estimates, digits=4)
#' @export
#' @family PanelSelect
#' @references
#' Bailey, M., & Peng, J. (2025). A Random Effects Model of Non-Ignorable Nonresponse in Panel Survey Data. Available at SSRN <https://www.ssrn.com/abstract=5475626>
probitRE_probitRE = function(probit1, probit2, id.name, data=NULL, par=NULL, method='BFGS', rho_off=FALSE, tau_off=FALSE, H=10, init=c('zero', 'unif', 'norm', 'default')[4], rho.init=0, tau.init=0, use.optim=FALSE, verbose=0){
  # 1.1 Sort data by id
  data = data.frame(data)
  data_original = copy(data)
  ord = order(data[, id.name])
  data = data[ord, ]
  group = c(0,cumsum(table(as.integer(factor(data[, id.name])))))

  # 1.1 parse outcome formula
  mf = model.frame(probit2, data=data, na.action=NULL, drop.unused.levels=TRUE)
  y = model.response(mf, "numeric")
  x = model.matrix(attr(mf, "terms"), data=mf)

  # 1.2 parse selection formula
  mf2 = model.frame(probit1, data=data, na.action=NULL, drop.unused.levels=TRUE)
  d = model.response(mf2, "numeric")
  w = model.matrix(attr(mf2, "terms"), data=mf2)

  # 1.3 Initialize parameters
  est_probit1 = glm(probit1, data=data, family=binomial(link="probit"))
  par_probit1 = coef(summary(est_probit1))[,1]
  est_probit2 = glm(probit2, data=data, family=binomial(link="probit"))
  par_probit2 = coef(summary(est_probit2))[,1]
  names(par_probit1) = paste0('probit1.', names(par_probit1))
  names(par_probit2) = paste0('probit2.', names(par_probit2))
  par_probit1[is.na(par_probit1)] = 0
  par_probit2[is.na(par_probit2)] = 0
  # convert bounded parameters to unboudned, rho to 2*atanh(rho)=ln((1+rho)/(1-rho))
  par = c(par_probit1, par_probit2, log_delta=0, log_lambda=0, artanh_rho=2*atanh(rho.init), artanh_tau=2*atanh(tau.init))
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
    res = optim(par=par, fn=LL_probitRE_probitRE, gr=Gradient_probitRE_probitRE, method=method, control=list(fnscale=-1, trace=as.numeric(verbose>=1)), y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off, verbose=verbose, hessian = TRUE)
    res$iterations = res$counts['function']
    res$LL = res$value
  }else{
    res = maxLik(LL_probitRE_probitRE, grad=Gradient_probitRE_probitRE, start=par, y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off, method=method, verbose=verbose, printLevel=verbose+1)
    res$par = res$estimate
    res$LL = res$maximum
  }
  res$n_obs = length(d)

  # 3. Compile results
  # res = getVarSE(res, verbose=verbose)
  gvar = Gradient_probitRE_probitRE(res$par,y,d,x,w,group,H,rho_off=rho_off,tau_off=tau_off,verbose=verbose-1,variance=TRUE)
  if(use.optim) res$gradient = gvar$g # optim has no gradient
  res = getVarSE(res, gvar=gvar, verbose=verbose)

  # res$num_g = numericGradient(LL_probitRE_probitRE,res$par,y=y, d=d, x=x, w=w, group=group, H=H, rho_off=rho_off, tau_off=tau_off)
  # cat('-------Gradient difference------\n')
  # print(res$num_g - gvar$g)

  # convert unbounded estimates to bounded estimates
  trans_vars=c(delta='log_delta', lambda='log_lambda', rho='artanh_rho', tau='artanh_tau')
  trans_types=c('exp', 'exp', 'correlation', 'correlation')
  if(rho_off){
    trans_vars = trans_vars[-(length(trans_vars)-1)]
    trans_types = trans_types[-(length(trans_types)-1)]
  }
  if(tau_off){
    trans_vars = trans_vars[-length(trans_vars)]
    trans_types = trans_types[-length(trans_types)]
  }
  res = transCompile(res, trans_vars, trans_types)
  res$predict = predict_probitRE_probitRE(res$estimates[, 1], res$var, probit1, probit2, data_original)

  # Need to estimate probitRE and probitRE models to make the test meaningful
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
