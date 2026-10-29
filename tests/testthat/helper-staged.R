# Shared settings for staged GL accuracy tests.
#
# Accuracy checks must not fail because a slow host exhausts the production
# wall-clock budgets (60 s per panel and 300 s per call by default). This keeps
# every numerical setting and gives each active GL option list generous
# explicit budgets. Budget enforcement itself is tested with mocked clocks in
# test-quadrature-staged.R.
staged_test_budget = function(control, seconds = 1e4) {
  for (field in c('quadrature_options', 'observation_options'))
    if (length(control[[field]]))
      control[[field]][c('seconds_per_panel', 'seconds_total')] = seconds
  control
}
