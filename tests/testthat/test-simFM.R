test_that("", {
  
  seed <- 01102025
  no_of_obs <- 102
  no_of_vars <- 50
  no_of_factors <- 3
  trans_error_var_cov <- diag(1, no_of_factors)
  loading_matrix <- matrix(round(rnorm(no_of_vars * no_of_factors)), no_of_vars, no_of_factors)
  meas_error_mean <- rep(0, no_of_vars)
  meas_error_var_cov <- diag(1, no_of_vars)
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
  factor_lag_order <- 2
  simul_delay <- c(floor(rexp(no_of_vars, 1)))
  quarterfy <- TRUE
  quarterly_variable_ratio  <- 0.10
  corr <- TRUE
  beta_param <- 2
  set.seed(seed)
  burn_in <- 999
  starting_date <- "1970-01-01"
  rescale <- TRUE
  check_stationarity <- TRUE
  stationarity_check_threshold <- 1e-10
  
  # Generate pseudo-date data
  expect_silent(data_date <- simFM(no_of_obs = no_of_obs, no_of_vars = no_of_vars,
                                   no_of_factors = no_of_factors, loading_matrix = loading_matrix,
                                   meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                                   trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff,
                                   factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy,
                                   quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr,
                                   beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                                   rescale = rescale, check_stationarity = check_stationarity,
                                   stationarity_check_threshold = stationarity_check_threshold))
  
  # Generate numeric data
  expect_silent(data_mat <- simFM(no_of_obs = no_of_obs, no_of_vars = no_of_vars,
                                  no_of_factors = no_of_factors, loading_matrix = loading_matrix,
                                  meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                                  trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff,
                                  factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy,
                                  quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr,
                                  beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = NULL,
                                  rescale = rescale, check_stationarity = check_stationarity,
                                  stationarity_check_threshold = stationarity_check_threshold))
  
  # Basic checks
  expect_equal(t(coredata(data_date$data)), data_mat$data)
  expect_equal(t(coredata(data_date$factors)), data_mat$factors)
  expect_equal(data_date$trans_var_coeff, data_mat$trans_var_coeff)
  expect_equal(data_date$loading_matrix, data_mat$loading_matrix)
  expect_equal(t(coredata(data_date$meas_error)), data_mat$meas_error)
  expect_equal(data_date$meas_error_var_cov, data_mat$meas_error_var_cov)
  expect_equal(data_date$trans_error_var_cov, data_mat$trans_error_var_cov)
  expect_equal(data_date$frequency, data_mat$frequency)
  expect_equal(data_date$delay, data_mat$delay)
  
})
