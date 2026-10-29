## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PanelSelect)

## ----old-call, eval=FALSE-----------------------------------------------------
# # PanelSelect 1.0.1
# old = probitRE_linearRE(d ~ x + w, y ~ x + w, id.name = "id", data = dat)

## ----current-call, eval=FALSE-------------------------------------------------
# # PanelSelect 2.0.0
# fit = panel_select(
#   outcome = y ~ x + w,
#   selection = d ~ x + w,
#   data = dat,
#   id = "id",
#   family = "linear")
# summary(fit)

