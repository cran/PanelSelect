#' Sample Selection Models for Panel Data
#' @description  This package supports a series of panel sample selection models, where the first stage is a panel Probit model with individual random effects and the second stage can be a panel linear, Probit, Poisson, or Poisson log-normal model with individual random effects. Models for count outcome are imported from the PanelCount package. \cr\cr
#' @section Functions:
#' probitRE_linearRE: panel sample selection model with continuous outcome \cr \cr
#' probitRE_probitRE: panel sample selection model with binary outcome \cr \cr
#' probitRE_PoissonRE: panel sample selection model with count outcome \cr \cr
#' probitRE_PLNRE: panel sample selection model with count outcome \cr \cr
#' @name PanelSelect
#' @importFrom statmod gauss.quad
#' @importFrom stats binomial rnorm dnorm pnorm qnorm dpois dlogis plogis model.frame model.matrix model.response optim pchisq poisson runif lm glm coef sigma logLik
#' @importFrom maxLik maxLik numericGradient numericHessian
#' @importFrom MASS mvrnorm
#' @importFrom pbivnorm pbivnorm
#' @importFrom pbv dbvnorm
#' @importFrom PanelCount ProbitRE_PoissonRE ProbitRE_PLNRE
#' @rawNamespace import(data.table)
#' @importFrom Rcpp sourceCpp
#' @useDynLib PanelSelect
NULL
'_PACKAGE'
#> NULL
