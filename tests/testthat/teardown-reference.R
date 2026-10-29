# Restore the exact registration state, including previously absent methods.
for (reference_method in names(reference_registration$previous)) {
  previous_method = reference_registration$previous[[reference_method]]
  if (is.null(previous_method))
    rm(list = reference_method, envir = reference_registration$table)
  else assign(reference_method, previous_method, envir = reference_registration$table)
}
rm(reference_registration, reference_method, previous_method)
