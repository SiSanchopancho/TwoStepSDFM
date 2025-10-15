# TwoStepSDFM
A ``C++``-based ``R`` implementation of a two-step estimation procedure for a (linear Gaussian) Sparse Dynamic Factor Model (SDFM) as outlined in Franjic and Schweikert (2024).

## Introduction

The ``TwoStepSDFM`` package provides a fast implementation of the Kalman Filter and Smoother (hereinafter KFS, see Koopman and Durbin, 2000) to estimate factors in a mixed-frequency SDFM framework, explicitly accounting for cross-sectional correlation in the measurement error. The KFS is initialized using results from Sparse Principal Components Analysis (SPCA) by Zou and Hastie (2006) in a preliminary step. This approach generalizes the two-step estimator for approximate dynamic factor models by Giannone, Reichlin, and Small (2008) and Doz, Giannone, and Reichlin (2011). For more details see Franjic and Schweikert (2024).

## Features

- **Fast Model Simulation**: The ``simFM()`` function provides a flexible framework to simulate approximate DFMs.
- **Fast Model Estimation**: The ``twoStepSDFM()`` function provides a fast and convenient implementation of the two-step estimator outlined in Franjic and Schweikert (2024).
- **Fast Model Prediction**: The ``nowcast()`` function is a highly convenient prediction function that automatically takes care of many issues that arise with mixed frequency data and ragged edges.
- **Compatibility**: All functions take advantage of ``C++`` for enhanced speed.

## Prerequisites

- **Rcpp**: A package for integrating `C++` code into `R` [@Eddelbuettel2011Rcpp]. [Rcpp CRAN repository](https://cran.r-project.org/web/packages/Rcpp/index.html)
- **RcppEigen**: A package for integrating the `Eigen` linear algebra library into `R` [@Bates2013EcppEigen]. [RcppEigen CRAN repository](https://cran.r-project.org/web/packages/RcppEigen/index.html)
- **GCC compiler** (version 5.0 or later) [GCC Website](https://gcc.gnu.org/).

## Installation

### Compile from scratch

``Rcpp`` and ``RcppEigen`` can be downloaded from CRAN or directly installed from within `R` by calling ``install.packages("...")``.

To install the package itself, a short `R` script is provided (see `PackageBuilder.R`). The package currently only compiles with the ``g++``/``gcc`` compiler. All tests have been performed and created using the ``C++14`` standard.

## Usage

### ``simFM()``: Simulating a dynamic factor model

#### Parameters

- ``no_of_observations`` Integer number of observations.
- ``no_of_variables`` Integer number of Variables.
- ``no_of_factors`` Integer number of factors.
- ``loading_matrix`` Numeric (no_of_variables x no_of_factors) loading matrix.
- ``meas_error_mean`` Numeric vector/matrix of the means of the measurement errors.
- ``meas_error_var_cov`` Numeric (no_of_factors x no_of_factors) varianceNumeric (no_of_variables x no_of_factors)-covariance matrix of the measurement errors.
- ``trans_error_var_cov`` Numeric (no_of_variables x no_of_variables) variance-covariance matrix of the transition errors.
- ``trans_var_coeff`` Either a list of length max_factor_lag_order with each entry a numeric (no_of_factors x no_of_factors) VAR coefficient matrix or a matrix of dimensions (no_of_factors x(no_of_factors * max_factor_lag_order)) holding the VAR coefficients of the factor VAR process in each (no_of_factors x no_of_factors) block.
- ``factor_lag_order`` Integer order of the VAR process in the state equation.
- ``delay`` Integer vector of delays imposed onto the end of the data in months (ragged edges)
- ``quarterfy`` Logical, whether or not some of the data should be aggregated to quarterly representations.
- ``quarterly_variable_ratio`` Ratio of variables ought to be quarterfied.
- ``corr`` Logical, whether or not the measurement error should be randomly correlated inside the function using a random correlation matrix with off-diagonal elements governed by a beta-distribution.
- ``beta_param`` Parameter of the beta-distribution governing the off-diagonal elements of the variance-covariance matrix of the measurement error.
- ``seed`` 32-bit unsigned integer seed for all random processes inside the function.
- ``burn_in`` Integer burn-in period of the simulated data ought to be discarded at the beginning of the data.
- ``rescale`` Logical, whether or not the variance of the measurement error should be rescaled by the common component to equalise the signal-to-noise ratio.
- ``starting_date`` A date type object indicating the start of the dataset. If NULL (default), the function returns matrices with observations along the second dimension (i.e., time in columns). If specified, the function treats the data as a time series, aligning it accordingly.
- ``check_stationarity`` Logical, whether or not the stationarity properties of the factor VAR process should be checked.
- ``stationarity_check_threshold`` Threshold of the stationarity check for when to deem an eigenvalue negative beyond the numerical error.
- ``parallel`` Logical, make use of Eigen internal parallel matrix operations

#### Example

```R
  set.seed(02102025)
  no_of_observations <- 200
  no_of_variables <- 150
  no_of_factors <- 3
  trans_error_var_cov <- diag(1, no_of_factors)
  loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors)
  meas_error_mean <- rep(0, no_of_variables)
  meas_error_var_cov <- diag(1, no_of_variables)
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
  factor_lag_order <- 2 # Order of the factor VAR process
  simul_delay <- c(6, 3, 6, 0, 3, rep(0, 45), floor(rexp(100, 1)))
  quarterfy <- TRUE
  quarterly_variable_ratio  <- 1/3
  corr <- TRUE
  beta_param <- 2 
  seed <- 01102025 
  set.seed(seed)
  burn_in <- 999 
  starting_date <- "1970-01-01"
  rescale <- TRUE
  check_stationarity <- TRUE
  stationarity_check_threshold <- 1e-10
  
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
              no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
              meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
              factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
              quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
              beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
              rescale = rescale, check_stationarity = check_stationarity, 
              stationarity_check_threshold = stationarity_check_threshold)
              
  print(FM)
  plot(FM)
```

### ``twoStepSDFM()``: Estimating the model parameters and higher-frequency factors

#### Parameters

- ``data`` Numeric (no_of_variables x no_of_observations) matrix of data or zoo/xts object.
- ``delay`` Integer vector of variable delays.
- ``selected`` Integer vector of the number of selected variables for each factor.
- ``no_of_factors`` Integer number of factors.
- ``max_factor_lag_order`` Integer max P of the VAR(P) process of the factors.
- ``decorr_errors`` Logical, whether or not the errors should be decorrelated.
- ``lag_estim_criterion`` Information criterion used for the estimation of the factor VAR order ("BIC", "AIC", "HIC").
- ``ridge_penalty`` Ridge penalty.
- ``lasso_penalty`` Numeric vector, lasso penalties for each factor (set to NaN if not intended as stopping criterion).
- ``max_iterations`` Integer maximum number of iterations.
- ``max_no_steps`` Integer number of max_no_steps (used for LARS-EN as an alternative).
- ``comp_null`` Computational zero.
- ``check_rank`` Logical, whether or not the rank of the variance-covariance matrix should be checked.
- ``conv_crit`` Conversion criterion for the SPCA algorithm.
- ``conv_threshold`` Conversion criterion for the coordinate descent algorithm.
- ``log`` Logical, whether or not output should be printed along the algorithm.
- ``parallel`` Logical, whether or not to use Eigen's internal parallel matrix operations.
- ``fcast_horizon`` Integer forecasting horizon for factor forecasts

#### Example

```R
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
  
  FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
              no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
              meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
              trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
              factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
              quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
              beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
              rescale = rescale, check_stationarity = check_stationarity, 
              stationarity_check_threshold = stationarity_check_threshold)
  
  frequency <- FM$frequency
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
  
  fit <- twoStepSDFM(data = FM$data, delay = delay, selected = selected, 
                     no_of_factors = no_of_factors,  max_factor_lag_order  = max_factor_lag_order,
                     decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                     ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                     max_iterations = max_iterations, max_no_steps = max_no_steps, 
                     comp_null = comp_null,  check_rank = check_rank,  conv_crit = conv_crit, 
                     conv_threshold = conv_threshold, log = log, parallel = parallel,
                     fcast_horizon = fcast_horizon)
                     
  print(fit)
  plot(fit)
```

### ``nowcast()``: Predict quarterly data using the two-step SDFM procedure of Franjic & Scheikert (2024).

#### Usage

Most of the parameters of ``nowcast()`` are directly parsed to ``twoStepSDFM()`` and work thus accordingly. The majour differences, however, lie in the parameter ``data``, ``variables_of_interest``, ``max_fcast_horizon``, ``frequency``, ``max_ar_lag_order``, and ``max_predictor_lag_order``.

``data`` expects a zoo/xts object of mixed frequency data. At least one variable in ``data`` must be of quarterly frequency. All quarterly variables should be stored such that their actual realisations are stored in the last months of each respective quarter. The inter-quarter month should be filled with either the realisation or some other non-NA numerical value. ``NA``s are fine at the end of the panel as long as they align with the delays provided in ``delay``. The function will check for ``NA``s that lie outside of the ragged edges.

``variables_of_interest`` is the vector indicating thos quarterly variables ought to be predicted. Note that currently only quarterly target variables are supported. In general, multiple quarterly series can be predicted at once.

``max_fcast_horizon`` indicates the maximum number of forecasts outside of the final observation of the panel should be computed. Note, the minimum forecasting horizon will be inferred internally for each target variable. Here, the minimum forecasting horizon will be set such that all ragged edges due to publication delay will be predicted. For example: Say the target variable is delayed by six months. ``max_fcast_horizon`` is set to two. The function will then automatically compute a one-step back backcast, a nowcast, a one-step ahead forecast, and a two step-ahead forecast. If the second variable of interest is not published with delay, only the one- and two-step ahead forecast will be computed.

#### Forecasting Method

The ``nowcast()`` function employs a forecast averaging forecasting method. In the first step, the factors are computed using the ``twoStepSDFM()`` function. Next, the computed factors are aggregated according to Mariano and Murasawa (2003) into a quarterly representation. They are then treated as quarterly predictors along the potential additional provided observed quarterly predictors. 

The function is generally able to compute predictions for multiple target variables at once. For each target variable and horizon, a prediction from an ARDL model is computed for each quarterly predictor and each of the aggregated factors. These predictions are then averaged (with equal weights) to a single prediction for each horizon.

#### Parameters

- ``data`` Zoo/xts object.
- ``variables_of_interest`` Integer vector indicating the index of all target variables
- ``max_fcast_horizon`` Maximum forecasting horizon
- ``delay`` Integer vector of predictor delays.
- ``selected`` Integer vector of the number of selected variables for each factor.
- ``frequency`` Integer vector of frequencies of the variables in the data set (currently supported: "12" for monthly and "4" for quarterly data)
- ``no_of_factors`` Integer number of factors.
- ``max_factor_lag_order`` Integer max P of the VAR(P) process of the factors.
- ``decorr_errors`` Logical, whether or not the errors should be decorrelated.
- ``lag_estim_criterion`` Information criterion used for the estimation of the factor VAR order ("BIC", "AIC", "HIC").
- ``ridge_penalty`` Ridge penalty.
- ``lasso_penalty`` Numeric vector, lasso penalties for each factor (set to NaN if not intended as stopping criterion).
- ``max_iterations`` Integer maximum number of iterations.
- ``max_no_steps`` Integer number of max_no_steps (used for LARS-EN as an alternative).
- ``comp_null`` Computational zero.
- ``check_rank`` Logical, whether or not the rank of the variance-covariance matrix should be checked.
- ``conv_crit`` Conversion criterion for the SPCA algorithm.
- ``conv_threshold`` Conversion criterion for the coordinate descent algorithm.
- ``log`` Logical, whether or not output should be printed along the algorithm.
- ``parallel`` Logical, whether or not to use Eigen's internal parallel matrix operations.
- ``max_ar_lag_order`` Integer maximum number of lags of the target variable ought to be included in the nowcasting step
- ``max_predictor_lag_order`` Integer maximum number of lags of the predictors ought to be included in the nowcasting step

#### Example

```R
  set.seed(02102025)
  no_of_observations <- 200
  no_of_variables <- 150
  no_of_factors <- 3
  trans_error_var_cov <- diag(1, no_of_factors)
  loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors)
  meas_error_mean <- rep(0, no_of_variables)
  meas_error_var_cov <- diag(1, no_of_variables)
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
  factor_lag_order <- 2 # Order of the factor VAR process
  simul_delay <- c(6, 3, 6, 0, 3, rep(0, 45), floor(rexp(100, 1)))
  quarterfy <- TRUE
  quarterly_variable_ratio  <- 1/3
  corr <- TRUE 
  beta_param <- 2 
  seed <- 01102025 
  set.seed(seed)
  burn_in <- 999 
  starting_date <- "1970-01-01"
  rescale <- TRUE
  check_stationarity <- TRUE
  stationarity_check_threshold <- 1e-10
  
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
  
  print(results)
  plot(results)
  
```

### ``crossVal()``: Validate the model hyperparameters using a random hyper-parameter-candidate-scheme (Bergstra and Bengio, 2021) via BIC (Despois and Doz, 2023) or time series cross-validation (Hyndman and Athanasopoulos,
2018).

#### Parameters

- ``data`` ``Zoo``/``xts`` object.
- ``variable_of_interest`` Integer value indicating the index of the single target variable
- ``fcast_horizon`` Integer value indicating the target horizon
- ``delay`` Integer vector of predictor delays.
- ``frequency`` Integer vector of frequencies of the variables in the data set (currently supported: ``12`` for monthly and ``4`` for quarterly data)
- ``no_of_factors`` Integer number of factors.
- ``seed`` 32-bit positive integer for drawing the random hyper-parameter candidates
- ``min_ridge_penalty`` Lower bound for the sampled ridge penalty coefficient
- ``max_ridge_penalty`` Upper bound for the sampled ridge penalty coefficient
- ``cv_repititions`` Integer number of predictions ought to be computed for each candidate set
- ``cv_size`` Integer number of candidate sets
- ``lasso_penalty_type`` Character indicating the lasso penalty type. ``"selected"`` uses the number of non-zero elements of the loading matrix. ``"penalty"`` uses the lasso size constraint directly. ``"steps"`` uses the number of steps.
- ``min_max_penalty`` Vector of size tow, where the first element indicates the lower and the second element indicates the upper bound of the lasso penalty equivalent. If lasso_penalty_type is set to ``"selected"`` or ``"steps"``, both elements must be strictly positive integers.
- ``max_factor_lag_order`` Integer max P of the VAR(P) process of the factors.
- ``decorr_errors`` Logical, whether or not the errors should be decorrelated.
- ``lag_estim_criterion`` Information criterion used for the estimation of the factor VAR order (``"BIC"``, ``"AIC"``, ``"HIC"``).
- ``ridge_penalty`` Ridge penalty.
- ``lasso_penalty`` Numeric vector, lasso penalties for each factor (set to ``NaN`` if not intended as stopping criterion).
- ``max_iterations`` Integer maximum number of iterations.
- ``max_no_steps`` Integer number of max_no_steps (used for LARS-EN as an alternative).
- ``comp_null`` Computational zero.
- ``check_rank`` Logical, whether or not the rank of the variance-covariance matrix should be checked.
- ``conv_crit`` Conversion criterion for the SPCA algorithm.
- ``conv_threshold`` Conversion criterion for the coordinate descent algorithm.
- ``parallel`` Logical, whether or not to run the cross-validation loop in parallel.
- ``max_ar_lag_order`` Integer maximum number of lags of the target variable ought to be included in the nowcasting step
- ``max_predictor_lag_order`` Integer maximum number of lags of the predictors ought to be included in the nowcasting step

#### Example

```{R}
  # Simulate a DGP using simFM
  set.seed(02102025)
  no_of_observations <- 102 + 3
  no_of_variables <- 50
  no_of_factors <- 2
  trans_error_var_cov <- diag(1, no_of_factors)
  loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors)
  meas_error_mean <- rep(0, no_of_variables)
  meas_error_var_cov <- diag(1, no_of_variables)
  trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
  factor_lag_order <- 2
  simul_delay <- c(3, floor(rexp(no_of_variables - 1, 1)))
  quarterfy <- TRUE
  quarterly_variable_ratio  <- 1/no_of_variables
  corr <- TRUE
  beta_param <- 2 
  seed <- 01102025 
  set.seed(seed)
  burn_in <- 999
  starting_date <- "1970-01-01"
  rescale <- TRUE
  check_stationarity <- TRUE
  stationarity_check_threshold <- 1e-10
  
  # Draw the FM object
  is_FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
                 no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
                 meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
                 trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
                 factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
                 quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
                 beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
                 rescale = rescale, check_stationarity = check_stationarity, 
                 stationarity_check_threshold = stationarity_check_threshold)
  
  data <- is_FM$data
  variable_of_interest <- 1
  fcast_horizon <- 0
  delay <- simul_delay
  frequency <- is_FM$frequency
  no_of_factors <- 3
  seed <- 09102025
  min_ridge_penalty <- 0.01
  max_ridge_penalty <- 1
  lasso_penalty_type <- "selected"
  min_max_penalty <- c(10, no_of_variables - 1)
  # lasso_penalty_type <- "steps"
  # min_max_penalty <- c(10, 50 * (no_of_variables - 1))
  # lasso_penalty_type <- "penalty"
  # min_max_penalty <- c(0.0001, 10)
  cv_repititions <- 3
  cv_size <- 100
  max_factor_lag_order = 10
  decorr_errors = TRUE
  lag_estim_criterion = "BIC"
  ridge_penalty = 1e-6
  lasso_penalty = NULL
  max_iterations = 100
  max_no_steps = NULL
  comp_null = 1e-15
  check_rank = FALSE
  conv_crit = 1e-4
  conv_threshold = 1e-4
  log = FALSE
  max_ar_lag_order = 5
  max_predictor_lag_order = 5
  
  cv <- crossVal(data = data, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                 delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                 seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                 cv_repititions = cv_repititions, cv_size = cv_size, lasso_penalty_type = lasso_penalty_type,
                 min_max_penalty = min_max_penalty, max_factor_lag_order = max_factor_lag_order,
                 decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                 ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, max_iterations = max_iterations,
                 max_no_steps = max_no_steps, comp_null = comp_null, check_rank = check_rank,
                 conv_crit = conv_crit, conv_threshold = conv_threshold, parallel = FALSE,
                 max_ar_lag_order = max_ar_lag_order, max_predictor_lag_order = max_predictor_lag_order)
  
  cv_parallel <- crossVal(data = data, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
                          delay = delay, frequency = frequency, no_of_factors = no_of_factors,
                          seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
                          cv_repititions = cv_repititions, cv_size = cv_size, lasso_penalty_type = lasso_penalty_type,
                          min_max_penalty = min_max_penalty, max_factor_lag_order = max_factor_lag_order,
                          decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                          ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, max_iterations = max_iterations,
                          max_no_steps = max_no_steps, comp_null = comp_null, check_rank = check_rank,
                          conv_crit = conv_crit, conv_threshold = conv_threshold, parallel = TRUE,
                          max_ar_lag_order = max_ar_lag_order, max_predictor_lag_order = max_predictor_lag_order)
```

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

© 2024-2025 Domenic Franjic

This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

This work is partially based on the LARS-EN and SPCA algorithms found in @Zou2020elasticnet.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

**To Contribute:**

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes with descriptive messages.
4. Push to your fork and submit a pull request.

## Support

If you have any questions or need assistance, please open an issue on the GitHub repository or contact us via email.

## Contact

- **Name**: Domenic Franjic
- **Institution**: University of Hohenheim
- **Department**: Econometrics and Statistics, Core Facility Hohenheim
- **E-Mail**: franjic@uni-hohenheim.de

### References

- Bates, Douglas, and Dirk Eddelbuettel. 2013. “Fast and Elegant Numerical Linear Algebra Using the RcppEigen Package.” Journal of Statistical Software 52 (5): 1–24.
- Bergstra, J. and Bengio, Y. (2012). "Random search for hyper-parameter optimization." Journal of Machine Learning Research, 13(2).
- Despois, T. and Doz, C. (2023). "Identifying and interpreting the factors in factor models via sparsity: Different approaches." Journal of Applied Econometrics, 38(4):533–555.
- Doz, Catherine, Domenico Giannone, and Lucrezia Reichlin. 2011. “A two-step estimator for large approximate dynamic factor models based on Kalman filtering.” Journal of Econometrics 164 (1): 188–205.
- Eddelbuettel, Dirk, and Romain François. 2011. “Rcpp: Seamless R and C++ Integration.” Journal of Statistical Software 40 (8): 1–18.
- Franjic, Domenic, and Karsten Schweikert. 2024. “Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model.” Available at SSRN.
- Giannone, Domenico, Lucrezia Reichlin, and David Small. 2008. “Nowcasting: The Real-Time Informational Content of Macroeconomic Data.” Journal of Monetary Economics 55 (4): 665–76.
- Guennebaud, Gaël, Benoît Jacob, et al. 2010. “Eigen V3.” http://eigen.tuxfamily.org.
- Hyndman, R. J. and Athanasopoulos, G. 2018. "Forecasting: principles and practice". OTexts Melbourne, 3 edition.
- Koopman, Siem Jan, and James Durbin. 2000. “Fast Filtering and Smoothing for Multivariate State Space Models.” Journal of Time Series Analysis 21 (3): 281–96.
- Mariano, Roberto S., and Yasutomo Murasawa. 2003. “A New Coincident Index of Business Cycles Based on Monthly and Quarterly Series.” Journal of Applied Econometrics 18 (4): 427–43.
- Zou, Hui, and Trevor Hastie. 2020. "Elasticnet: Elastic-Net for Sparse Estimation and Sparse PCA."" https://CRAN.R-project.org/package=elasticnet.
- Zou, Hui, Trevor Hastie, and Robert Tibshirani. 2006. “Sparse Principal Component Analysis.” Journal of Computational and Graphical Statistics 15 (2): 265–86.

