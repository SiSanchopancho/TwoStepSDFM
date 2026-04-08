test_that("", {
  
  # Load data
  data("mixed_freq_factor_model")
  data_zoo <- scale(mixed_freq_factor_model$data)
  no_of_factors <- dim(mixed_freq_factor_model$factors)[2]
  no_of_mtly_vars <- sum(mixed_freq_factor_model$frequency == 12)
  variable_of_interest <- 1
  fcast_horizon <- 1
  delay <- mixed_freq_factor_model$delay
  frequency <- mixed_freq_factor_model$frequency
  seed <- 09102025
  min_ridge_penalty <- 0.01
  max_ridge_penalty <- 1
  lasso_penalty_type <- "selected"
  min_max_penalty <- c(10, no_of_mtly_vars - 1)
  cv_repetitions <- 3
  cv_size <- 10
  parallel = TRUE
  no_of_cores <- 2
  max_ar_lag_order = 5
  max_predictor_lag_order = 5
  
  # Set seed to check whether cross-val changes with the global RNG stat
  set.seed(seed)
  old_seed <- .Random.seed
  
  # CV results in series
  expect_no_error(cv_series <- crossVal(data = data_zoo, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                                        delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                                        seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                                        cv_repetitions = cv_repetitions, cv_size = cv_size, lasso_penalty_type = lasso_penalty_type,
                                        min_max_penalty = min_max_penalty, verbose = FALSE))
  
  # CV results in parallel
  expect_no_error(cv_parallel <- crossVal(data = data_zoo, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                                          delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                                          seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                                          cv_repetitions = cv_repetitions, cv_size = cv_size, lasso_penalty_type = lasso_penalty_type,
                                          min_max_penalty = min_max_penalty, verbose = FALSE, parallel = TRUE, no_of_cores = no_of_cores))
  
  # Basic checks (set call null as they are different by construction)
  cv_series$call <- NULL
  cv_parallel$call <- NULL
  expect_equal(cv_series, cv_parallel)
  
  
  # Test selecting according to the max number of steps
  expect_no_error(cv_series_steps <- crossVal(data = data_zoo, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                                              delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                                              seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                                              cv_repetitions = cv_repetitions, cv_size = cv_size, lasso_penalty_type = "steps",
                                              min_max_penalty = c(1, 500), verbose = FALSE))
  expect_no_error(cv_parallel_steps <- crossVal(data = data_zoo, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                                                delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                                                seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                                                cv_repetitions = cv_repetitions, cv_size = cv_size, lasso_penalty_type = "steps",
                                                min_max_penalty = c(1, 500), parallel = TRUE, 
                                                no_of_cores = no_of_cores, verbose = FALSE))
  cv_series_steps$call <- NULL
  cv_parallel_steps$call <- NULL
  expect_equal(cv_series_steps, cv_parallel_steps)
  
  # Test selecting according to the lasso penalty
  expect_no_error(cv_series_lasso <- crossVal(data = data_zoo, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                                              delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                                              seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                                              cv_repetitions = cv_repetitions, cv_size = cv_size, lasso_penalty_type = "penalty",
                                              min_max_penalty = c(0.0001, 10), verbose = FALSE))
  expect_no_error(cv_parallel_lasso <- crossVal(data = data_zoo, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                                                delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                                                seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                                                cv_repetitions = cv_repetitions, cv_size = cv_size, lasso_penalty_type = "penalty",
                                                min_max_penalty = c(0.0001, 10), parallel = TRUE, 
                                                no_of_cores = no_of_cores, verbose = FALSE))
  cv_series_lasso$call <- NULL
  cv_parallel_lasso$call <- NULL
  expect_equal(cv_series_lasso, cv_parallel_lasso)
  
  # Check whether cross-val changes with the global RNG stat
  current_seed <- .Random.seed
  expect_equal(old_seed, current_seed)
  
})
