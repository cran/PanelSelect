# Independent comparison code stays outside the production namespace.
# reference/math.R supplies independent base-R likelihood and score calculations.
reference_environment = new.env(parent = environment())
for (reference_file in c('derivatives.R', 'math.R', 'mixture.R', 'numerics.R', 'pln-inner.R',
  'prediction.R', 'quadrature.R', 'staged.R'))
  sys.source(testthat::test_path('reference', reference_file), envir = reference_environment)
reference_methods = mget(c('Ops.ps_jet', 'Math.ps_jet'), envir = reference_environment)
for (reference_name in setdiff(ls(reference_environment), names(reference_methods)))
  assign(reference_name, get(reference_name, envir = reference_environment), envir = environment())
rm(reference_environment, reference_file, reference_name)
