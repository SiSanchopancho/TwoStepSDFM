test_that("", {
  
  # Simulate a DGP using simFM
  set.seed(02102025)
  no_of_observations <- 200 # Number of observations
  no_of_variables <- 150 # Number of variabes
  no_of_factors <- 3 # Number of factors
  trans_error_var_cov <- diag(1, no_of_factors) # Variance-covariance matrix of the transition errors
  loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors) # Factor loadings matrix
  meas_error_mean <- rep(0, no_of_variables) # Mean of the measurement error
  meas_error_var_cov <- diag(1, no_of_variables) # Variance-covariance matrix of the measurement error
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors)) # Factor VAR coefficient matrix
  factor_lag_order <- 2 # Order of the factor VAR process
  simul_delay <- c(6, 3, 6, 0, 3, rep(0, 45), floor(rexp(100, 1)))
  quarterfy <- TRUE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
  quarterly_variable_ratio  <- 1/3
  corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
  beta_param <- 2 # Beta parameter governing the degree of correlation of the measurement error
  seed <- 01102025 # Seed
  set.seed(seed)
  burn_in <- 999 # Burn-in period
  starting_date <- "1970-01-01"
  rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 
  check_stationarity <- TRUE
  stationarity_check_threshold <- 1e-10
  
  # Simulate data
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
              no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
              meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
              factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
              quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
              beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
              rescale = rescale, check_stationarity = check_stationarity, 
              stationarity_check_threshold = stationarity_check_threshold)
  
  data <- FM$data
  variables_of_interest <- 1:5
  max_fcast_horizon <- 4
  delay <- simul_delay
  selected <- c(round(no_of_variables * 0.5), round(no_of_variables * 0.25))
  frequency <- c(rep(4, 50), rep(12, 100))
  no_of_factors <- 2 
  max_factor_lag_order  <- 10
  decorr_errors <- TRUE 
  lag_estim_criterion  <- "BIC"
  ridge_penalty <- 1e-06
  lasso_penalty <- NULL
  max_iterations <- 1000
  max_no_steps <- NULL
  comp_null <- 1e-15
  check_rank <- FALSE
  conv_crit <- 1e-04 
  conv_threshold <- 1e-04
  log <- FALSE
  parallel <- FALSE
  max_ar_lag_order <- 5
  max_predictor_lag_order  <- 5
  
  # Minimal working example
  expect_silent(nowcast(data = data, variables_of_interest = variables_of_interest, 
                        max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                        frequency = frequency, no_of_factors = no_of_factors))
  
  # Test single factor
  expect_silent(nowcast(data = data, variables_of_interest = variables_of_interest, 
                        max_fcast_horizon = max_fcast_horizon, delay = delay, selected = no_of_factors,
                        frequency = frequency, no_of_factors = 1))
  
  # Test single variable of interest
  expect_silent(nowcast(data = data, variables_of_interest = 1, 
                        max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                        frequency = frequency, no_of_factors = no_of_factors))
  
  # Test single quarterly predictor
  expect_silent(nowcast(data = data[, c(1, 51:150)], variables_of_interest = 1, 
                        max_fcast_horizon = max_fcast_horizon, delay = delay[c(1, 51:150)], selected = selected,
                        frequency = frequency[c(1, 51:150)], no_of_factors = no_of_factors))
  
  # Test single quarterly predictor but several targets
  expect_silent(nowcast(data = data[, c(1:3, 51:150)], variables_of_interest = 1:3, 
                        max_fcast_horizon = max_fcast_horizon, delay = delay[c(1:3, 51:150)], selected = selected,
                        frequency = frequency[c(1:3, 51:150)], no_of_factors = no_of_factors))
  
  # One factor, one target, no additional predictors
  expect_silent(nowcast(data = data[, c(1, 51:150)], variables_of_interest = 1, 
                        max_fcast_horizon = max_fcast_horizon, delay = delay[c(1, 51:150)], selected = no_of_factors,
                        frequency = frequency[c(1, 51:150)], no_of_factors = 1))
  
  
  results <- nowcast(data = data, variables_of_interest = variables_of_interest, 
                     max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                     frequency = frequency, no_of_factors = no_of_factors, 
                     max_factor_lag_order  = max_factor_lag_order,  
                     decorr_errors = decorr_errors, 
                     lag_estim_criterion  = lag_estim_criterion,
                     ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                     max_iterations = max_iterations, max_no_steps = max_no_steps, 
                     comp_null = comp_null, check_rank = check_rank,  
                     conv_crit = conv_crit, conv_threshold = conv_threshold, 
                     log = log, parallel = parallel, max_ar_lag_order = max_ar_lag_order,
                     max_predictor_lag_order = max_predictor_lag_order
  )
  
  # Basic checks
  expect_equal(sum(results$`SDFM Fit`$loading_matrix_estimate == 0), 87)
  expect_type(results, "list")
  expect_true(all(c("Forecasts", "Single Predictor Forecasts", "SDFM Fit")
                  %in% names(results)))
  
  results_determinism_check <- nowcast(data = data, variables_of_interest = variables_of_interest, 
                                       max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                                       frequency = frequency, no_of_factors = no_of_factors, 
                                       max_factor_lag_order  = max_factor_lag_order,  
                                       decorr_errors = decorr_errors, 
                                       lag_estim_criterion  = lag_estim_criterion,
                                       ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                                       max_iterations = max_iterations, max_no_steps = max_no_steps, 
                                       comp_null = comp_null, check_rank = check_rank,  
                                       conv_crit = conv_crit, conv_threshold = conv_threshold, 
                                       log = log, parallel = parallel, max_ar_lag_order = max_ar_lag_order,
                                       max_predictor_lag_order = max_predictor_lag_order
  )
  
  expect_equal(results_determinism_check, results)
  
  # Misuse (Note: Most parameter misuse is already checked in test-twoStepSDFM.R)
  expect_error(nowcast(data = data, variables_of_interest = -1, 
                       max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                       frequency = frequency, no_of_factors = no_of_factors, 
                       max_factor_lag_order  = max_factor_lag_order,  
                       decorr_errors = decorr_errors, 
                       lag_estim_criterion  = lag_estim_criterion,
                       ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                       max_iterations = max_iterations, max_no_steps = max_no_steps, 
                       comp_null = comp_null, check_rank = check_rank,  
                       conv_crit = conv_crit, conv_threshold = conv_threshold, 
                       log = log, parallel = parallel, max_ar_lag_order = max_ar_lag_order,
                       max_predictor_lag_order = max_predictor_lag_order
                       ))
  
  expect_error(nowcast(data = data, variables_of_interest = variables_of_interest, 
                       max_fcast_horizon = -1, delay = delay, selected = selected,
                       frequency = frequency, no_of_factors = no_of_factors, 
                       max_factor_lag_order  = max_factor_lag_order,  
                       decorr_errors = decorr_errors, 
                       lag_estim_criterion  = lag_estim_criterion,
                       ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                       max_iterations = max_iterations, max_no_steps = max_no_steps, 
                       comp_null = comp_null, check_rank = check_rank,  
                       conv_crit = conv_crit, conv_threshold = conv_threshold, 
                       log = log, parallel = parallel, max_ar_lag_order = max_ar_lag_order,
                       max_predictor_lag_order = max_predictor_lag_order
  ))
  
  expect_error(nowcast(data = data, variables_of_interest = variables_of_interest, 
                       max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                       frequency = frequency, no_of_factors = no_of_factors, 
                       max_factor_lag_order  = max_factor_lag_order,  
                       decorr_errors = decorr_errors, 
                       lag_estim_criterion  = lag_estim_criterion,
                       ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                       max_iterations = max_iterations, max_no_steps = max_no_steps, 
                       comp_null = comp_null, check_rank = check_rank,  
                       conv_crit = conv_crit, conv_threshold = conv_threshold, 
                       log = log, parallel = parallel, max_ar_lag_order = -1,
                       max_predictor_lag_order = max_predictor_lag_order
  ))
  
  expect_error(nowcast(data = data, variables_of_interest = variables_of_interest, 
                       max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                       frequency = frequency, no_of_factors = no_of_factors, 
                       max_factor_lag_order  = max_factor_lag_order,  
                       decorr_errors = decorr_errors, 
                       lag_estim_criterion  = lag_estim_criterion,
                       ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                       max_iterations = max_iterations, max_no_steps = max_no_steps, 
                       comp_null = comp_null, check_rank = check_rank,  
                       conv_crit = conv_crit, conv_threshold = conv_threshold, 
                       log = log, parallel = parallel, max_ar_lag_order = max_ar_lag_order,
                       max_predictor_lag_order = -1
  ))
  
  
})