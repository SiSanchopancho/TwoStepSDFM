test_that("Simulaiton function simFM works", {
  
  # Draw from an approximate mixed frequency factor model
  T <- 100 # Number of observations
  N <- 5 # Number of variabes
  R <- 2 # Number of factors
  Sigma_epsilon <- diag(1, R) # Variance-covariance matrix of the transition errors
  Lambda <- matrix(rnorm(N * R), N, R) # Factor loadings matrix
  mu_xi <- rep(0, N) # Mean of the measurement error
  Sigma_xi <- diag(1, N) # Variance-covariance matrix of the measurement error
  Phi <- cbind(diag(0.5, R), -diag(0.25, R)) # Factor VAR coefficient matrix
  P <- 2 # Order of the factor VAR process
  quarterfy <- TRUE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
  corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
  beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
  m <- 0.03 # Ratio of monthly predictors ought to be quarterfied
  seed <- 16022024 # Seed
  set.seed(seed)
  burn_in <- 999 # Burn-in period
  rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 
  
  # Draw data
  FM <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
                 Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = quarterfy, m = m,
                 corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale)
  
  expect_type(FM, "list")
  expect_true(all(c("factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error_var_cov",
                    "transition_error_var_cov", "measurement_error", "data", "frequency") %in% names(FM)))
  expect_equal(dim(FM$data), c(N, T))
  expect_equal(dim(FM$factors), c(R, T))
  
  FM_determinism_check <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
              Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = quarterfy, m = m,
              corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale)
  expect_equal(FM$data, FM_determinism_check$data, tolerance = 1e-10)
  
  
})

test_that("Estimation function twoStepSDFM works", {
  
  
  # Simulate a DGP using simFM
  T <- 100 # Number of observations
  N <- 5 # Number of variabes
  R <- 2 # Number of factors
  Sigma_epsilon <- diag(1, R) # Variance-covariance matrix of the transition errors
  Lambda <- matrix(rnorm(N * R), N, R) # Factor loadings matrix
  mu_xi <- rep(0, N) # Mean of the measurement error
  Sigma_xi <- diag(1, N) # Variance-covariance matrix of the measurement error
  Phi <- cbind(diag(0.5, R), -diag(0.25, R)) # Factor VAR coefficient matrix
  P <- 2 # Order of the factor VAR process
  quarterfy <- FALSE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
  m <- 0
  corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
  beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
  seed <- 16022024 # Seed
  set.seed(seed)
  burn_in <- 999 # Burn-in period
  rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 
  
  # Draw the FM object
  FM <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
              Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = quarterfy, m = m,
              corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale,
              check_stationarity = TRUE, stationarity_check_threshold = 1e-10)
  
  # Fitting a sparse model with l2 regularisation and non-orthogonal measurement errors
  selected <- c(round(N * 0.8), round(N * 0.5))
  delay <- rep(0, N)
  fit <- twoStepSDFM(FM$data, delay, selected, R, l2 = 1e-4, parallel = FALSE)
  
  expect_type(fit, "list")
  expect_true(all(c("loading_matrix_estimate", "filtered_state_variance", "factor_estimate",
                      "smoothed_state_variance", "error_var_cov_cholesky_factor",
                      "factor_var_lag_order") %in% names(fit)))
  expect_equal(dim(fit$factor_estimate), c(R, T + 1))
  
  # Draw the FM object
  fit_determinism_chec <- twoStepSDFM(FM$data, delay, selected, R, l2 = 1e-4, , parallel = FALSE)
  expect_equal(fit$factor_estimate, fit_determinism_chec$factor_estimate, tolerance = 1e-4)
  expect_equal(fit$loading_matrix_estimate, fit_determinism_chec$loading_matrix_estimate, tolerance = 1e-4)
  
  
})
