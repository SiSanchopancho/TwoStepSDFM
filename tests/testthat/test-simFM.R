test_that("simFM() works", {
  
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
  
  # Minimal example
  expect_silent(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                      no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                      meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                      trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                      factor_lag_order = factor_lag_order))
  
  # Test single factor
  expect_silent(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                      no_of_factors = 1, loading_matrix = matrix(rnorm(no_of_variables), no_of_variables, 1), 
                      meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                      trans_error_var_cov = matrix(1, 1, 1), trans_var_coeff = matrix(c(0.5, -0.2), 1, 2), 
                      factor_lag_order = factor_lag_order))
  
  # Test correct usage with pseudo dates
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
              no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
              meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
              factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
              quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
              beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
              rescale = rescale, check_stationarity = check_stationarity, 
              stationarity_check_threshold = stationarity_check_threshold)
  
  expect_type(FM, "list")
  expect_true(all(c("data", "factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error",
                    "measurement_error_var_cov", "transition_error_var_cov", "transition_error_var_cov",
                    "frequency", "delay") %in% names(FM)))
  expect_equal(dim(FM$data), c(no_of_observations, no_of_variables))
  expect_equal(dim(FM$factors), c(no_of_observations, no_of_factors))
  FM_determinism_check <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                                no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                                meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                                trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                                factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                                quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                                beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                                rescale = rescale, check_stationarity = check_stationarity, 
                                stationarity_check_threshold = stationarity_check_threshold)
  expect_equal(FM$data, FM_determinism_check$data, tolerance = 1e-10)
  expect_true(is.zoo(FM$data))
  expect_true(is.zoo(FM$factors))
  expect_true(is.zoo(FM$measurement_error))
  
  # Test correct usage without pseudo dates
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
              no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
              meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
              factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
              quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
              beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = NULL,
              rescale = rescale, check_stationarity = check_stationarity, 
              stationarity_check_threshold = stationarity_check_threshold)
  
  expect_type(FM, "list")
  expect_true(all(c("data", "factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error",
                    "measurement_error_var_cov", "transition_error_var_cov", "transition_error_var_cov",
                    "frequency", "delay") %in% names(FM)))
  expect_equal(dim(FM$data), c(no_of_variables, no_of_observations))
  expect_equal(dim(FM$factors), c(no_of_factors, no_of_observations))
  FM_determinism_check <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                                no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                                meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                                trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                                factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                                quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                                beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = NULL,
                                rescale = rescale, check_stationarity = check_stationarity, 
                                stationarity_check_threshold = stationarity_check_threshold)
  expect_equal(FM$data, FM_determinism_check$data, tolerance = 1e-10)
  expect_true(!is.zoo(FM$data))
  expect_true(!is.zoo(FM$factors))
  expect_true(!is.zoo(FM$measurement_error))
  
  # Test mishandling
  expect_error(simFM(no_of_observations = -1, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = -1, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = -1, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = t(loading_matrix), 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = 1:2, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = -1,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = -1, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = -1, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = -1, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = 1:2, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = 5, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = -1, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = 2^63, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = -1, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = "This is not a date",
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = -8, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = FLASE, 
                     stationarity_check_threshold = stationarity_check_threshold))
  
  expect_error(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                     no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                     meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                     trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                     factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                     quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                     beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                     rescale = rescale, check_stationarity = check_stationarity, 
                     stationarity_check_threshold = -1))
  
  # Check warnings
  unit_root_matrix <- diag(no_of_factors)
  expect_warning(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                       no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                       meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                       trans_error_var_cov = trans_error_var_cov, trans_var_coeff = unit_root_matrix, 
                       factor_lag_order = 1, delay = simul_delay, quarterfy = quarterfy, 
                       quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                       beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                       rescale = rescale, check_stationarity = TRUE, 
                       stationarity_check_threshold = 1e-10))
  
  non_stationary <- cbind(diag(no_of_factors), diag(no_of_factors))
  expect_warning(simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                       no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                       meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                       trans_error_var_cov = trans_error_var_cov, trans_var_coeff = non_stationary, 
                       factor_lag_order = 2, delay = simul_delay, quarterfy = quarterfy, 
                       quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                       beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                       rescale = rescale, check_stationarity = TRUE, 
                       stationarity_check_threshold = 1e-10))
  
})
