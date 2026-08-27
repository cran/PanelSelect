#' @inherit PanelCount::ProbitRE_PLNRE title description details return references
#' @inheritParams PanelCount::ProbitRE_PLNRE
#' @note This function is imported from the *PanelCount* package (see \link[PanelCount]{ProbitRE_PLNRE} for details).
#' @examples
#' \donttest{
#' library(MASS)
#' library(PanelSelect)
#' set.seed(1)
#' N = 400
#' periods = 5
#' rho = 0.5
#' tau = 0
#'
#' id = rep(1:N, each=periods)
#' time = rep(1:periods, N)
#' x = rnorm(N*periods)
#' w = rnorm(N*periods)
#'
#' # correlated random effects at the individual level
#' r = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
#' r1 = rep(r[,1], each=periods)
#' r2 = rep(r[,2], each=periods)
#'
#' # correlated error terms at the individual-time level
#' e = mvrnorm(N*periods, mu=c(0,0), Sigma=matrix(c(1,tau,tau,1), nrow=2))
#' e1 = e[,1]
#' e2 = e[,2]
#'
#' # selection
#' z = as.numeric(1+x+w+r1+e1>0)
#' # outcome
#' y = rpois(N*periods, exp(-1+x+r2+e2))
#' y[z==0] = NA
#' dt = data.frame(id,time,x,w,z,y)
#'
#' # As N increases, the parameter estimates will be more accurate
#' m = probitRE_PLNRE(z~x+w, y~x, data=dt, id.name='id', verbose=-1)
#' print(m$estimates, digits=4)
#' }
#' @export
#' @family PanelSelect
#' @importFrom PanelCount ProbitRE_PLNRE
probitRE_PLNRE <- ProbitRE_PLNRE

#' @inherit PanelCount::ProbitRE_PoissonRE title description details return examples references
#' @inheritParams PanelCount::ProbitRE_PoissonRE
#' @note This function is imported from the *PanelCount* package (see \link[PanelCount]{ProbitRE_PoissonRE} for details).
#' @examples
#' library(MASS)
#' library(PanelSelect)
#' set.seed(1)
#' N = 500
#' periods = 5
#' rho = 0.5
#'
#' id = rep(1:N, each=periods)
#' time = rep(1:periods, N)
#' x = rnorm(N*periods)
#' w = rnorm(N*periods)
#'
#' # correlated random effects at the individual level
#' r = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
#' r1 = rep(r[,1], each=periods)
#' r2 = rep(r[,2], each=periods)
#'
#' e = rnorm(N*periods)
#'
#' # selection
#' z = as.numeric(1+x+w+r1+e>0)
#' # outcome
#' y = rpois(N*periods, exp(-1+x+r2))
#' y[z==0] = NA
#' dt = data.frame(id,time,x,w,z,y)
#'
#' # As N increases, the parameter estimates will be more accurate
#' m = probitRE_PoissonRE(z~x+w, y~x, data=dt, id.name='id', verbose=-1)
#' print(m$estimates, digits=4)
#' @export
#' @family PanelSelect
#' @importFrom PanelCount ProbitRE_PoissonRE
probitRE_PoissonRE <- ProbitRE_PoissonRE
