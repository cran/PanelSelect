# Keep registration test-only, using APIs available at the testthat minimum.
reference_registration = local({
  table = get('.__S3MethodsTable__.', envir = asNamespace('base'))
  methods = names(reference_methods)
  previous = mget(methods, envir = table, inherits = FALSE, ifnotfound = list(NULL))
  for (method in methods)
    registerS3method(sub('\\.ps_jet$', '', method), 'ps_jet', reference_methods[[method]])
  list(table = table, previous = previous)
})
