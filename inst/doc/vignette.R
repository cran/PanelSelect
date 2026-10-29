## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PanelSelect)

## ----data---------------------------------------------------------------------
set.seed(100)
N = 500L
period = 5L
n = N*period
id = rep(1:N, each = period)
u = rnorm(N)
v = -0.5*u + sqrt(0.75)*rnorm(N)
eps = rnorm(n)
e = 0.5*eps + sqrt(0.75)*rnorm(n)
w = rnorm(n)
x = rnorm(n)
d = as.integer(x + w + u[id] + eps > 0)
y = x + w + v[id] + e
y[d == 0] = NA
dat = data.frame(id = id, time = rep(seq_len(period), N), x = x, w = w,
                 d = d, y = y)

## ----selection-fit------------------------------------------------------------
fit = panel_select(
  outcome = y ~ x + w,
  selection = d ~ x + w,
  data = dat,
  id = "id",
  family = "linear",
  control = panel_control(quadrature = "gh", nodes = 10L))

## ----numerical-checks---------------------------------------------------------
check_stationarity(fit)[c("optimizer_converged", "Q_OPG")]

## ----quadrature-check---------------------------------------------------------
check_quadrature(fit, nodes = 20L)

## ----refine-integration, eval=FALSE-------------------------------------------
# fit_more = update(fit, control = modifyList(fit$control, list(nodes = 20L)))
# coef(fit_more) - coef(fit)

