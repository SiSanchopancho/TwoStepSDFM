#' ## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib TwoStepSDFM
## usethis namespace: end
NULL

# SPDX-License-Identifier: GPL-3.0-or-later
#
#  Copyright \u00A9 2024 Domenic Franjic
#
#  This file is part of TwoStepSDFM.
#
#  TwoStepSDFM is free software: you can redistribute
#  it and/or modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  TwoStepSDFM is distributed in the hope that it
#  will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.


#' @name twoStepSDFM
#' @title Run Sparse DFM KFS
#'
#' This function wraps the C++ implementation of the Sparse DFM KFS to
#' provide default parameters and potentially pre-process inputs or
#' post-process outputs.
#'
#' @param data Numeric "no_of_variables x no_of_observations" matrix of data.
#' @param delay Integer vector of predictor delays.
#' @param selected Integer vector of the number of selected variables for each factor.
#' @param no_of_factors Integer, number of factors.
#' @param max_factor_lag_order Integer, max P of the VAR(P) process of the factors.
#' @param decorr_errors Logical, whether or not the errors should be decorrelated.
#' @param lag_estim_criterion Character string, information criterion used for the estimation of the factor VAR order ("BIC", "AIC", "HIC").
#' @param ridge_penalty Numeric, ridge penalty.
#' @param lasso_penalty Numeric vector, lasso penalties for each factor (set to NaN if not intended as stopping criterion).
#' @param max_iterations Integer, maximum number of iterations.
#' @param max_no_steps Integer, number of max_no_steps (used for LARS-EN as an alternative).
#' @param comp_null Numeric, computational zero.
#' @param check_rank Logical, whether or not the rank of the variance-covariance matrix should be checked.
#' @param conv_crit Numeric, conversion criterion for the SPCA algorithm.
#' @param conv_threshold Numeric, conversion criterion for the coordinate descent algorithm.
#' @param log Logical, whether or not output should be printed along the algorithm.
#' @param parallel Make use of Eigen internal parallel matrix operations
#' @return 
#' The `twoStepSDFM` function returns an object `fit` which, contains the following elements:
#' \describe{
#'   \item{\code{Lambda_hat}}{A matrix of estimated loadings for each factor on each observed variable.}
#'   \item{\code{F}}{The factor estimates.}
#'   \item{\code{Pt}}{The estimated precision matrix of the factors.}
#'   \item{\code{Wt}}{The estimated variance-covariance matrix of the residuals/errors.}
#'   \item{\code{C}}{A matrix representing any transformations applied to the factors for error decorrelation, if \code{decorr_errors} is TRUE.}
#'   \item{\code{max_factor_lag_order}}{The order of the VAR(max_factor_lag_order) model used in the factor model process.}
#' }
#' @export
twoStepSDFM <- function(data,
                        delay,
                        selected,
                        no_of_factors,
                        max_factor_lag_order = 10,
                        decorr_errors = TRUE,
                        lag_estim_criterion = "BIC",
                        ridge_penalty = 1e-6,
                        lasso_penalty = NaN,
                        max_iterations = 1000,
                        max_no_steps = NaN,
                        comp_null = 1e-15,
                        check_rank = FALSE,
                        conv_crit = 1e-4,
                        conv_threshold = 1e-4,
                        log = FALSE,
                        parallel = FALSE) {
  
  no_of_variables <- dim(data)[1]
  no_of_observations <- dim(data)[1]
  
  # Misshandling of the loadings matrix
  
  data_r <- try(as.matrix(data), silent = TRUE)
  if (inherits(data_r, "try-error")) {
    stop(paste0("data must be a matrix or convertable to a matrix"))
  }
  
  data_r <- t(data_r)
  
  # Misshandling of selected
  
  if (length(selected) == 1) {
    if (is.nan(selected) || is.null(selected) || is.infinite(-selected)) {
      selected <- as.integer(rep(dim(data_r)[2], no_of_factors))
    }
  }else{
    selected <- as.integer(selected)
  }
  
  if (length(selected) != no_of_factors || any(selected < 0)) {
    stop(paste0("The vector containing the number of variables ought to be selected is either not of length no_of_factors = ", no_of_factors, " or has negative values. If you aim to disable selected as stopping criterion set it to NaN, -Inf, or Null."))
  }
  
  if (length(delay) == 1) {
    if (is.nan(delay) || is.null(delay) || is.infinite(-delay)) {
      delay <- as.integer(rep(0, no_of_variables))
    }
  }else{
    delay <- as.integer(delay)
  }
  
  if (length(delay) != no_of_variables || any(delay < 0)) {
    stop(paste0("The vector containing the monthly delays is either not of length no_of_variables = ", no_of_variables, " or has negative values. If you aim to provide no delay set it to NaN, -Inf, or Null."))
  }
  
  l1_disabled <- is.null(lasso_penalty) || (length(lasso_penalty) == 1 && (is.nan(lasso_penalty) || is.infinite(-lasso_penalty)))
  if (!l1_disabled && (length(lasso_penalty) != no_of_factors || any(lasso_penalty < 0))) {
    stop(paste0("The vector of lasso_penalty penalties is either not of length no_of_factors = ", no_of_factors, " or has negative values. If you aim to disable the lasso_penalty penalty as stopping criterin set it to NaN, -Inf, or Null."))
  }
  
  if (l1_disabled) {
    lasso_penalty <- rep(-2147483647L, no_of_factors)
  }
  
  if (!is.na(max_no_steps) && !is.null(max_no_steps)) {
    if (max_no_steps < 0 && !is.infinite(-max_no_steps)) {
      stop(paste0("The maximum number of max_no_steps used for the LARS-EN algorith must be non-negative. To disable the number of max_no_steps as stopping criterion, set \"max_no_steps\" to NaN, -Inf, or Null"))
    }
  }
  
  if (is.nan(max_no_steps) || is.null(max_no_steps) || is.infinite(-max_no_steps)) {
    max_no_steps <- -2147483647L # C++ INT_MIN
  }
  
  if (max_no_steps != floor(max_no_steps) && !is.na(max_no_steps)) {
    stop(paste0("The maximum number of max_no_steps used for the LARS-EN algorith must be an integer."))
  }
  max_no_steps <- as.integer(max_no_steps)
  
  # Disable the checking of the KFS converions
  
  KFS_conv_crit <- 1e+15
  
  # Call the C++ function
  
  result <- runSDFMKFS(
    X_in = data_r,
    delay = delay,
    selected = selected,
    R = as.integer(no_of_factors),
    order = as.integer(max_factor_lag_order),
    decorr_errors = decorr_errors,
    crit = lag_estim_criterion,
    l2 = ridge_penalty,
    l1 = lasso_penalty,
    max_iterations = as.integer(max_iterations),
    steps = max_no_steps,
    comp_null = comp_null,
    check_rank = check_rank,
    conv_crit = conv_crit,
    conv_threshold = conv_threshold,
    log = log,
    KFS_conv_crit = KFS_conv_crit,
    parallel = parallel
  )
  
  # Rename the results
  names(result) <- c("loading_matrix_estimate", "filtered_state_variance", "factor_estimate",
                     "smoothed_state_variance", "error_var_cov_cholesky_factor",
                     "factor_var_lag_order")
    
  # Return the result
  return(result)
}


#' @name simFM
#' @title Draw data from a dynamic factor model
#'
#' The function simulates an (approximate) Dynamic Factor Model
#'
#' @param no_of_observations Number of observations.
#' @param no_of_variables Number of Variables.
#' @param no_of_factors Number of factors.
#' @param loading_matrix Loadings matrix.
#' @param meas_error_mean Mean of the measurement errors.
#' @param meas_error_var_cov Variance-covariance matrix of the measurement errors.
#' @param trans_error_var_cov Variance-covariance matrix of the transition errors.
#' @param trans_var_coeff Either a list of length max_factor_lag_order or a matrix of dimensions Rx(no_of_factors * max_factor_lag_order) or (no_of_factors * max_factor_lag_order)xR holding the VAR coefficients of the factor VAR process.
#' @param factor_lag_order Order of the VAR process in the state equation.
#' @param quarterfy Whether or not some of the data should be aggregated to quarterly representations.
#' @param quarterly_variable_ratio Ratio of variables ought to be quarterfied.
#' @param corr Whether or not the measurement error should be randomly correlated inside the function using a random correlation matrix with off-diagonal elements governed by a beta-distribution.
#' @param beta_param Parameter of the beta-distribution governing the off-diagonal elements of the variance-covariance matrix of the measurement error.
#' @param seed Seed for all random processes inside the function.
#' @param burn_in Burn-in period of the simulated data ought to be discarded at the beginning of the data.
#' @param rescale Whether or not the variance of the measurement error should be rescaled by the common component to equalise the signal-to-noise ratio.
#' @param check_stationarity Whether or not the stationarity properties of the factor VAR process should be checked.
#' @param stationarity_check_threshold Threshold of the stationarity check for when to deem an eigenvalue negative beyond the numerical error.
#' @param parallel Make use of Eigen internal parallel matrix operations
#' @return Returns an `simFM` object containing the DFM parameters expressed in the factor VAR companion form.
#' @export
simFM <- function(no_of_observations, no_of_variables, no_of_factors, loading_matrix, meas_error_mean, meas_error_var_cov, trans_error_var_cov, trans_var_coeff, factor_lag_order, quarterfy, quarterly_variable_ratio = 0, corr = FALSE,
                  beta_param = Inf, seed = 20022024, burn_in = 1000, rescale = TRUE,
                  check_stationarity = FALSE, stationarity_check_threshold = 1e-5,
                  parallel = FALSE) {
  # Check for miss-handling
  
  if(!is.logical(quarterfy) && !(quarterfy %in% c(0, 1))){
    
    stop(paste0("quarterfy should be logical."))
    
  }
  
  if(!is.logical(rescale) && !(rescale %in% c(0, 1))){
    
    stop(paste0("rescale should be logical."))
    
  }
  
  # Misshandling of the State-Space Dimensions
  
  if (!is.numeric(no_of_variables) || !is.numeric(no_of_observations) || !is.numeric(no_of_factors) || !is.numeric(factor_lag_order)) {
    stop(paste0("Either the number of variables no_of_variables, the number of Observations no_of_observations, the number of Factors no_of_factors, or the factor VAR order factor_lag_order is not provided as numeric."))
  }
  
  if (no_of_variables <= 0 || no_of_observations <= 0 || no_of_factors <= 0 || factor_lag_order <= 0) {
    stop(paste0("Either the number of variables no_of_variables, the number of Observations no_of_observations, the number of Factors no_of_factors, or the factor VAR process order factor_lag_order is negative or zero."))
  }
  
  if (max(c(no_of_variables, no_of_observations)) <= no_of_factors) {
    warning(paste0("no_of_factors = ", no_of_factors, " >= ", max(no_of_variables, no_of_observations), " = max(no_of_variables,no_of_observations)! In general, no_of_factors << no_of_variables,no_of_observations."))
  }
  
  # Misshandling of the loadings matrix
  
  LambdaR <- try(as.matrix(loading_matrix), silent = TRUE)
  if (inherits(LambdaR, "try-error")) {
    stop(paste0("loading_matrix must be a matrix or convertable to a matrix"))
  }
  
  if (dim(LambdaR)[1] != no_of_variables || dim(LambdaR)[2] != no_of_factors) {
    stop(paste0("The loadings Matrix loading_matrix must be of dimensions NxR = ", no_of_variables, "x", no_of_factors, " but is ", dim(LambdaR)[1], "x", dim(LambdaR)[2], "."))
  }
  
  # Misshandling of the measurement error mean
  
  muR <- try(as.matrix(meas_error_mean), silent = TRUE)
  if (inherits(muR, "try-error")) {
    stop(paste0("meas_error_mean must be a vector or convertable to a matrix"))
  }
  
  if (dim(muR)[1] == 1 && dim(muR)[2] == no_of_variables) {
    muR <- t(muR)
  }
  
  if(!is.numeric(muR)){
    
    stop(paste0("meas_error_mean must be of type numeric"))
    
  }
  
  if (dim(muR)[1] != no_of_variables && dim(muR)[2] != 1) {
    stop(paste0("The mean of the measurement error meas_error_mean must be a vector of dimensions ", no_of_variables, "x1 or 1x", no_of_variables, " but is ", dim(muR)[1], "x", dim(muR)[2], "."))
  }
  
  # Misshandling of the measurement error variance-covariance matrix
  
  SigmaR <- try(as.matrix(meas_error_var_cov), silent = TRUE)
  if (inherits(SigmaR, "try-error")) {
    stop(paste0("meas_error_var_cov must be a matrix or convertable to a matrix"))
  }
  
  if(!is.numeric(SigmaR)){
    
    stop(paste0("meas_error_var_cov must be of type numeric"))
    
  }
  
  if (sum(dim(SigmaR) == c(no_of_variables, no_of_variables)) != 2) {
    stop(paste0("The variance-covariance matrix of the measurement error meas_error_var_cov must be a vector of dimensions ", no_of_variables, "x", no_of_variables, " but is ", dim(SigmaR)[1], "x", dim(SigmaR)[2], "."))
  }
  
  if (!isSymmetric(SigmaR)) {
    stop(paste0("The variance-covariance matrix of the measurement error must be symmetric."))
  }
  
  # Misshandling of the transition error variance-covariance matrix
  
  SR <- try(as.matrix(trans_error_var_cov), silent = TRUE)
  if (inherits(SR, "try-error")) {
    stop(paste0("trans_error_var_cov must be a matrix or convertable to a matrix"))
  }
  
  if (!all(dim(SR) == c(no_of_factors, no_of_factors))) {
    stop(paste0("The variance-covariance matrix of the transition error trans_error_var_cov must be a vector of dimensions ", no_of_factors, "x", no_of_factors, " but is ", dim(SigmaR)[1], "x", dim(SigmaR)[2], "."))
  }
  
  if (!isSymmetric(SR)) {
    stop(paste0("The variance-covariance matrix of the transition error must be symmetric."))
  }
  
  if (!isSymmetric(SR)) {
    stop(paste0("The variance-covariance matrix of the transition error must be symmetric."))
  }
  
  # Misshandling of the VAR coefficient matrix
  
  if (is.list(trans_var_coeff) && !is.data.frame(trans_var_coeff)) {
    # If trans_var_coeff is provided as list
    
    # Check whether each element in trans_var_coeff has the correct dimensions and whether there is the correct number of matrices provided
    
    ind <- c()
    for (i in 1:length(trans_var_coeff)) {
      s <- sum(dim(trans_var_coeff[[i]]) == c(no_of_factors, no_of_factors))
      if (s != 2) {
        ind[i] <- i
      }
    }
    
    if (!is.null(ind)) {
      stop("The VAR coefficient matrices in trans_var_coeff must be of dimensions ", no_of_factors, "x", no_of_factors, ", which is not the case for the matrices ", paste(ind, collapse = ", "), ".")
    }
    
    if (length(trans_var_coeff) != factor_lag_order) {
      stop("The number of VAR coefficient matrices in trans_var_coeff must equal factor_lag_order = ", factor_lag_order, " but is ", length(trans_var_coeff), ".")
    }
    
    # Store the list elements as matrix
    
    PhiR <- matrix(0, no_of_factors, no_of_factors * factor_lag_order)
    for (o in 1:factor_lag_order) {
      PhiR[1:no_of_factors, ((o - 1) * no_of_factors + 1):((o - 1) * no_of_factors + no_of_factors)] <- trans_var_coeff[[o]]
    }
  } else {
    # If trans_var_coeff is provided as matrix or equivalent
    
    PhiR <- try(as.matrix(trans_var_coeff), silent = TRUE)
    if (inherits(PhiR, "try-error")) {
      stop(paste0("PhiR must be a list, a matrix or convertable to a matrix"))
    }
    
    if (all(dim(PhiR) == c(factor_lag_order * no_of_factors, no_of_factors))) {
      PhiR <- t(PhiR)
    }
    
    if (!all(dim(PhiR) == c(no_of_factors, no_of_factors * factor_lag_order))) {
      stop(paste0("The VAR coefficient matrix of the factor process must be a vector of dimensions ", no_of_factors, "x", no_of_factors * factor_lag_order, " or ", no_of_factors * factor_lag_order, "x", no_of_factors, ", where each ", no_of_factors, "x", no_of_factors, " block represents a VAR coefficient of corresponding lag factor_lag_order. Here, the provided matrix is ", dim(PhiR)[1], "x", dim(PhiR)[2], "."))
    }
  }
  
  if(!is.numeric(PhiR)){
    
    stop("trans_var_coeff must be a numeric matrix or a list of numeric matrices")
    
  }
  
  # Check whether the process will result in a stationary process
  
  if (check_stationarity) {
    
    Comp <- matrix(0, no_of_factors * factor_lag_order, no_of_factors * factor_lag_order)
    
    Comp[1:no_of_factors, 1:(factor_lag_order * no_of_factors)] <- PhiR
    
    Comp[(no_of_factors + 1):(factor_lag_order * no_of_factors), 1:(no_of_factors * (factor_lag_order - 1))] <- diag(1, no_of_factors * (factor_lag_order - 1))
    
    # Compute eigenvalues
    e <- eigen(Comp)$values
    test <- any(abs(abs(e) - 1) < stationarity_check_threshold)
    
    if (test) {
      warning("At least one eigenvalue very close to one detected. Process might have random walk properties.\n")
    }
    
    test2 <- all(abs(e) < 1 + stationarity_check_threshold)
    
    if (!test2) {
      warning("At least one eigenvalue lies outside the complex unit-circle. Process is not stationary.")
    }
    
  }
  
  # Misshandling of the beta parameter
  
  if (beta_param <= 0) {
    stop("The value of the parameter of the beta distribution governing the cross-sectional correlation in the measurement error beta_param must be strictly positive. To disable the random cross-section correlation of the measurement error set corr = FALSE.")
  }
  
  if (corr && is.infinite(beta_param)) {
    stop("If the data should be randomly correlated internally beta_param should not be Inf. To disable the random cross-section correlation of the measurement error set corr = FALSE.")
  }
  
  # M;isshandling of the ratio of quarterly variables
  
  if (!is.numeric(quarterly_variable_ratio)) {
    stop("The ratio of variables ought to be quarterfied quarterly_variable_ratio must be numeric.")
  }
  
  if (quarterly_variable_ratio < 0) {
    stop("The ratio of variables ought to be quarterfied quarterly_variable_ratio must be non-negative.")
  }
  
  # Misshandling of the burn_in period
  
  if (burn_in < 0) {
    stop("The burn in period at the start of the simulated data burn_in must be non-negative.")
  }
  
  FM <- runStaticFM(
    T = no_of_observations, N = no_of_variables, S = SR, Lambda = LambdaR, mu_e = muR, Sigma_e = SigmaR, A = PhiR, order = factor_lag_order, quarterfy = quarterfy, corr = corr, beta_param = beta_param,
    m = quarterly_variable_ratio, seed = seed, R = no_of_factors, burn_in = burn_in, rescale = rescale, parallel = parallel
  )
  
  names(FM) <- c("factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error_var_cov",
                 "transition_error_var_cov", "measurement_error", "data", "frequency")
  
  return(FM)
}
