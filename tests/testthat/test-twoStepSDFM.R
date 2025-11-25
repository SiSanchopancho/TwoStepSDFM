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
  
  frequency <- FM_no_date$frequency
  delay <- simul_delay
  selected <- c(10, 8, 5)
  no_of_factors <- 3 
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
  fcast_horizon <- 10
  
  # Minimal example
  
  # Test single factor
  expect_silent(twoStepSDFM(data = FM_date$data, delay = delay, selected = no_of_factors,
                            no_of_factors = 1))
  
  # Parcing zoo object
  expect_silent(twoStepSDFM(data = FM_date$data, delay = delay, selected = selected,
                            no_of_factors = no_of_factors))
  
  # Parcing matrix object
  expect_silent(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected,
                            no_of_factors = no_of_factors))
  
  fit_sparse_date <- twoStepSDFM(data = FM_date$data, delay = delay, selected = selected, 
                                 no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                                 decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                                 ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                                 max_iterations = max_iterations, max_no_steps = max_no_steps, 
                                 comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                                 conv_threshold = conv_threshold, log = log, parallel = parallel,
                                 fcast_horizon = fcast_horizon)
  
  fit_sparse_no_date <- twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                                    no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                                    decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                                    ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                                    max_iterations = max_iterations, max_no_steps = max_no_steps, 
                                    comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                                    conv_threshold = conv_threshold, log = log, parallel = parallel,
                                    fcast_horizon = fcast_horizon)
  
  # Basic checks
  expect_equal(sum(fit_sparse_date$loading_matrix_estimate == 0), 7)
  expect_equal(sum(fit_sparse_no_date$loading_matrix_estimate == 0), 7)
  expect_type(fit_sparse_date, "list")
  expect_type(fit_sparse_no_date, "list")
  expect_true(all(c("data", "loading_matrix_estimate", "smoothed_factors", "smoothed_state_variance",
                    "factor_var_lag_order", "error_var_cov_cholesky_factor")
                  %in% names(fit_sparse_date)))
  expect_true(all(c("data", "loading_matrix_estimate", "smoothed_factors", "smoothed_state_variance",
                    "factor_var_lag_order", "error_var_cov_cholesky_factor")
                  %in% names(fit_sparse_no_date)))
  expect_equal(dim(fit_sparse_date$smoothed_factors), c(no_of_observations + fcast_horizon, no_of_factors))
  expect_equal(dim(fit_sparse_no_date$smoothed_factors), c(no_of_factors, no_of_observations + fcast_horizon))
  fit_sparse_coredata <- t(coredata(fit_sparse_date$smoothed_factors))
  expect_equal(fit_sparse_coredata, fit_sparse_no_date$smoothed_factors, tolerance = 1e-10)
  expect_equal(fit_sparse_date$loading_matrix_estimate, fit_sparse_no_date$loading_matrix_estimate, tolerance = 1e-10)
  expect_true(is.zoo(fit_sparse_date$smoothed_factors))
  expect_true(!is.zoo(fit_sparse_no_date$smoothed_factors))
  
  # test information criteria
  expect_silent(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                            no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                            decorr_errors = decorr_errors,  lag_estim_criterion  = "AIC", 
                            ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                            max_iterations = max_iterations, max_no_steps = max_no_steps, 
                            comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                            conv_threshold = conv_threshold, log = log, parallel = parallel,
                            fcast_horizon = fcast_horizon))
  expect_silent(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                            no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                            decorr_errors = decorr_errors,  lag_estim_criterion  = "HIC", 
                            ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                            max_iterations = max_iterations, max_no_steps = max_no_steps, 
                            comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                            conv_threshold = conv_threshold, log = log, parallel = parallel,
                            fcast_horizon = fcast_horizon))
  
  # Test misshandling
  expect_error(twoStepSDFM(data = matrix(NaN, no_of_factors, no_of_observations), delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = 1, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = 1, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = -1,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = -5,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = "decorr_errors",  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = "CIB", 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = -5, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = NaN, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = -1, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = -55, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = -1,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = "FLASE",  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = NaN, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = -8, log = log, parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = "No Please", parallel = parallel,
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = "Yes Please",
                           fcast_horizon = fcast_horizon))
  
  expect_error(twoStepSDFM(data = FM_no_date$data, delay = delay, selected = selected, 
                           no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                           decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                           ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                           max_iterations = max_iterations, max_no_steps = max_no_steps, 
                           comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                           conv_threshold = conv_threshold, log = log, parallel = parallel,
                           fcast_horizon = -1))
  
  
})

