test_that("Simulaiton function simFM works", {
  
  # Simulate a DGP using simFM
  no_of_observations <- 100 # Number of observations
  no_of_variables <- 10 # Number of variabes
  no_of_factors <- 2 # Number of factors
  trans_error_var_cov <- diag(1, no_of_factors) # Variance-covariance matrix of the transition errors
  loading_matrix <- matrix(rnorm(no_of_variables * no_of_factors), no_of_variables, no_of_factors) # Factor loadings matrix
  meas_error_mean <- rep(0, no_of_variables) # Mean of the measurement error
  meas_error_var_cov <- diag(1, no_of_variables) # Variance-covariance matrix of the measurement error
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors)) # Factor VAR coefficient matrix
  factor_lag_order <- 2 # Order of the factor VAR process
  quarterfy <- FALSE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
  quarterly_variable_ratio  <- 0
  corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
  beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
  seed <- 16022024 # Seed
  set.seed(seed)
  burn_in <- 999 # Burn-in period
  rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 
  
  # Draw the FM object
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, no_of_factors = no_of_factors, loading_matrix = loading_matrix, meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, factor_lag_order = factor_lag_order, quarterfy = quarterfy, quarterly_variable_ratio  = quarterly_variable_ratio ,
              corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale,
              check_stationarity = TRUE, stationarity_check_threshold = 1e-10)

  
  expect_type(FM, "list")
  expect_true(all(c("factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error_var_cov",
                    "transition_error_var_cov", "measurement_error", "data", "frequency") %in% names(FM)))
  expect_equal(dim(FM$data), c(no_of_variables, no_of_observations))
  expect_equal(dim(FM$factors), c(no_of_factors, no_of_observations))
  
  FM_determinism_check <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, no_of_factors = no_of_factors, loading_matrix = loading_matrix, meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                                trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, factor_lag_order = factor_lag_order, quarterfy = quarterfy, quarterly_variable_ratio  = quarterly_variable_ratio ,
                                corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale,
                                check_stationarity = TRUE, stationarity_check_threshold = 1e-10)
  expect_equal(FM$data, FM_determinism_check$data, tolerance = 1e-10)

})

test_that("Estimation function twoStepSDFM works", {
  
  # Simulate a DGP using simFM
  no_of_observations <- 100 # Number of observations
  no_of_variables <- 10 # Number of variabes
  no_of_factors <- 2 # Number of factors
  trans_error_var_cov <- diag(1, no_of_factors) # Variance-covariance matrix of the transition errors
  loading_matrix <- matrix(rnorm(no_of_variables * no_of_factors), no_of_variables, no_of_factors) # Factor loadings matrix
  meas_error_mean <- rep(0, no_of_variables) # Mean of the measurement error
  meas_error_var_cov <- diag(1, no_of_variables) # Variance-covariance matrix of the measurement error
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors)) # Factor VAR coefficient matrix
  factor_lag_order <- 2 # Order of the factor VAR process
  quarterfy <- FALSE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
  quarterly_variable_ratio  <- 0
  corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
  beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
  seed <- 16022024 # Seed
  set.seed(seed)
  burn_in <- 999 # Burn-in period
  rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 
  
  # Draw the FM object
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, no_of_factors = no_of_factors, loading_matrix = loading_matrix, meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, factor_lag_order = factor_lag_order, quarterfy = quarterfy, quarterly_variable_ratio  = quarterly_variable_ratio,
              corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale,
              check_stationarity = TRUE, stationarity_check_threshold = 1e-10)

  selected <- c(round(no_of_variables * 0.8), round(no_of_variables * 0.5))
  delay <- rep(0, no_of_variables)
  fit_sparse <- twoStepSDFM(data = FM$data, delay = delay, selected = selected, no_of_factors = no_of_factors, 
                            max_factor_lag_order  = 10, decorr_errors = TRUE, 
                            lag_estim_criterion  = "BIC", ridge_penalty = 1e-06, 
                            lasso_penalty = NaN, max_iterations = 1000, max_no_steps = NaN, 
                            comp_null = 1e-15,  check_rank = FALSE,  conv_crit = 1e-04, 
                            conv_threshold = 1e-04, log = FALSE, parallel = FALSE)
  
  expect_type(fit_sparse, "list")
  expect_true(all(c("loading_matrix_estimate", "filtered_state_variance", "factor_estimate",
                      "smoothed_state_variance", "error_var_cov_cholesky_factor",
                      "factor_var_lag_order") %in% names(fit_sparse)))
  expect_equal(dim(fit_sparse$factor_estimate), c(no_of_factors * fit_sparse$factor_var_lag_order, no_of_observations))
  
  # Draw the FM object
  fit_determinism_chec <- twoStepSDFM(data = FM$data, delay = delay, selected = selected, no_of_factors = no_of_factors, 
                                      max_factor_lag_order  = 10, decorr_errors = TRUE, 
                                      lag_estim_criterion  = "BIC", ridge_penalty = 1e-06, 
                                      lasso_penalty = NaN, max_iterations = 1000, max_no_steps = NaN, 
                                      comp_null = 1e-15,  check_rank = FALSE,  conv_crit = 1e-04, 
                                      conv_threshold = 1e-04, log = FALSE, parallel = FALSE)
  expect_equal(fit_sparse$factor_estimate, fit_determinism_chec$factor_estimate, tolerance = 1e-4)
  expect_equal(fit_sparse$loading_matrix_estimate, fit_determinism_chec$loading_matrix_estimate, tolerance = 1e-4)
  
  
})
