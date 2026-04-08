test_that("", {
  
  # Load data
  data("mixed_freq_factor_model")
  data_zoo <- scale(mixed_freq_factor_model$data)
  no_of_factors <- dim(mixed_freq_factor_model$factors)[2]
  no_of_mtly_vars <- sum(mixed_freq_factor_model$frequency == 12)
  delay <- mixed_freq_factor_model$delay
  frequency <- mixed_freq_factor_model$frequency
  max_ar_lag_order <- 5
  max_predictor_lag_order <- 5
  variables_of_interest <- 1:2
  max_fcast_horizon <- 4
  selected <- rep(round(no_of_mtly_vars * 0.5), no_of_factors)
  
  # Parsing zoo object
  expect_silent(nowcast(data = data_zoo, variables_of_interest = variables_of_interest, 
                        max_fcast_horizon = 0, delay = delay, selected = selected,
                        frequency = frequency, no_of_factors = no_of_factors))
  
})
