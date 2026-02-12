test_that("", {
  
  # Simulate a DGP using simFM
  set.seed(02102025)
  no_of_observations <- 200
  no_of_variables <- 10
  no_of_factors <- 3
  trans_error_var_cov <- diag(1, no_of_factors)
  loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors)
  meas_error_mean <- rep(0, no_of_variables)
  meas_error_var_cov <- diag(1, no_of_variables)
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
  factor_lag_order <- 2 
  simul_delay <- c(floor(rexp(10, 1)))
  quarterfy <- FALSE
  quarterly_variable_ratio  <- 0
  corr <- TRUE
  beta_param <- 2
  seed <- 01102025
  set.seed(seed)
  burn_in <- 999
  starting_date <- "1970-01-01"
  rescale <- TRUE
  check_stationarity <- TRUE
  stationarity_check_threshold <- 1e-10
  
  FM_date <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                   no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                   meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                   trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                   factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                   quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                   beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                   rescale = rescale, check_stationarity = check_stationarity, 
                   stationarity_check_threshold = stationarity_check_threshold)
  
  FM_no_date <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                      no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                      meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                      trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                      factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                      quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                      beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = NULL,
                      rescale = rescale, check_stationarity = check_stationarity, 
                      stationarity_check_threshold = stationarity_check_threshold)
  
  
  # Minimal example
  
  # Test single factor
  expect_no_error(noOfFactors(data = FM_date$data, min_no_factors = 1, max_no_factors = 2))
  
  # Parsing zoo object
  expect_no_error(noOfFactors(data = FM_date$data, min_no_factors = 1, max_no_factors = 5))
  
  # Parsing matrix object
  expect_no_error(noOfFactors(data = FM_no_date$data, min_no_factors = 1, max_no_factors = 5))
  
  # Mishandling
  expect_error(noOfFactors(data = FM_no_date$data, min_no_factors = -1, max_no_factors = 5))
  expect_error(noOfFactors(data = FM_no_date$data, min_no_factors = 5, max_no_factors = 4))
  expect_error(noOfFactors(data = FM_no_date$data, min_no_factors = 1, max_no_factors = 5, 
                           confidence_threshold = 1.1))
  expect_error(noOfFactors(data = FM_no_date$data, min_no_factors = 1, max_no_factors = 5, 
                           confidence_threshold = -0.1))
  
})
