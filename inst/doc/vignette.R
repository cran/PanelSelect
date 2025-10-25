## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Table, echo=FALSE--------------------------------------------------------
tb = data.frame(
    Model = c('probitRE_linearRE', 'probitRE_probitRE', 'probitRE_PoissonRE', 'probitRE_PLNRE'),
    `First Stage` = c(rep('probitRE', 4)),
    `Second Stage` = c('linearRE', 'ProbitRE', 'PoissonRE', 'PLN_RE'),
    `Outcome Variable` = c('linear','binary', 'count', 'count'),
    Package = c('PanelSelect', 'PanelSelect', 'PanelCount', 'PanelCount')
)
knitr::kable(tb, col.names = gsub("[.]", " ", names(tb)), caption='**Table 1. Panel Sample Selection Models**')

## -----------------------------------------------------------------------------
library(PanelSelect)
library(MASS)
N = 500
period = 5
obs = N*period
rho = -0.5
tau = 0.5
set.seed(100)

re = mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), nrow=2))
u = rep(re[,1], each=period)
v = rep(re[,2], each=period)
e = mvrnorm(obs, mu=c(0,0), Sigma=matrix(c(1,tau,tau,1), nrow=2))
e1 = e[,1]
e2 = e[,2]

t = rep(1:period, N)
id = rep(1:N, each=period)
w = rnorm(obs)
z = rnorm(obs)
x = rnorm(obs)
d = as.numeric(x + w + u + e1 > 0)
y = x + w + v + e2
y[d==0] = NA
dt = data.frame(id, t, y, x, w, z, d)

# As N increases, the parameter estimates will be more accurate
m = probitRE_linearRE(d~x+w, y~x+w, 'id', dt, H=10, verbose=-1)
print(m$estimates, digits=4)

