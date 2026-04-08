#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @import zoo
#' @import xts
#' @import lubridate
#' @import ggplot2
#' @import stats
#' @import utils
#' @import patchwork
NULL

# SPDX-License-Identifier: GPL-3.0-or-later
#
#  Copyright (C) 2024-2026 Domenic Franjic
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

#' @name kalmanFilterSmoother
#' @title Univariate Representation of the Multivariate Kalman Filter and 
#' Smoother
#' @description
#' Filter and smooth the latent states/factors of a linear Gaussian state-space 
#' model, with measurement equation
#' \deqn{
#'   \bm{x}_t = \bm{\Lambda} \bm{f}_{t} + \bm{\xi}_t,\quad \bm{\xi}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{\xi}),
#' }
#' and transition equation
#' \deqn{
#'   \bm{f}_t = \sum_{p=0}^P\bm{\Phi}_p \bm{f}_{t-p} + \bm{\epsilon}_t,\quad \bm{\epsilon}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{f}).
#' }
#' for t = 1, \dots, T. For filtering and smoothing, the univariate 
#' representation of the multivariate Kalman Filter and Smoother is implemented 
#' according to \insertRef{koopman2000fast}{TwoStepSDFM}.
#' 
#' @param data Numeric (no_of_vars \eqn{\times}{x} no_of_obs) matrix of data or 
#' zoo/xts object sampled at the same frequency.
#' @param delay Integer vector of variable delays.
#' @param no_of_factors Integer number of factors.
#' @param loading_matrix Numeric (no_of_vars \eqn{\times}{x} no_of_factors) 
#' loading matrix.
#' @param meas_error_var_cov Numeric (no_of_factors \eqn{\times}{x} 
#' no_of_factors) variance-covariance matrix of the measurement errors.
#' @param trans_error_var_cov Numeric (no_of_vars \eqn{\times}{x} no_of_vars) 
#' variance-covariance matrix of the transition errors.
#' @param trans_var_coeff Either a list of length max_factor_lag_order with each
#' entry a numeric (no_of_factors \eqn{\times}{x} no_of_factors) VAR 
#' coefficient matrix or a matrix of dimensions (no_of_factors x(no_of_factors 
#' * max_factor_lag_order)) holding the VAR coefficients of the factor VAR 
#' process in each (no_of_factors \eqn{\times}{x} no_of_factors) block.
#' @param factor_lag_order Integer order of the VAR process in the state 
#' equation.
#' @param fcast_horizon Integer number of additional Filter predictions into the 
#' future.
#' @param decorr_errors Logical, whether or not the errors should be 
#' decorrelated (should be `TRUE` if `meas_error_var_cov` is not diagonal).
#' @param comp_null Computational zero.
#' @param parallel Logical, whether or not to use Eigen's internal parallel 
#' matrix operations.
#' @param jitter Numerical jitter for stability of internal solver algorithms. 
#' The jitter is added to the diagonal entries of the variance-covariance matrix 
#' of the measurement errors.
#' 
#' @details
#' To implement the univariate representation of the Kalman Filter and Smoother, 
#' the measurement error term has to be cross-sectionally uncorrelated. If 
#'`meas_error_var_cov` is not diagonal, one should set `decorr_errors = TRUE` so
#' that the data can be decorrelated internally prior to filtering and smoothing.
#'  
#' When decorrelating, the function first adds `jitter` to the diagonal elements 
#' of `meas_error_var_cov` and then tries to compute the Cholesky factor via 
#' Eigen's standard LLT decomposition \insertCite{eigenweb}{TwoStepSDFM}. If the
#' initial decorrelation fails, it silently switches to Eigen's more robust, 
#' but slower, LDLT decomposition with pivoting 
#' \insertCite{eigenweb}{TwoStepSDFM}. If this also fails, it is likely that 
#' `meas_error_var_cov` is not well-behaved. The analysis should be repeated 
#' with a larger `jitter` or a more robust variance-covariance matrix 
#' (estimator). The success of the internal Cholesky decomposition is reported 
#' by `llt_success_code`.
#' 
#' @return
#' An object of class `KFSFit` with components:
#' \describe{
#'   \item{data}{Original data matrix.}
#'   \item{smoothed_factors}{Object containing the smoothed factor estimates. 
#'   The object inherits its class from `data`: If `data` is provided as `zoo`, 
#'   `smoothed_factors` will be a `zoo` object. If `data` is provided as 
#'   `matrix`, `smoothed_factors` will be a (`no_of_factors` \eqn{\times}{x} 
#'   `no_of_obs`) matrix.}
#'   \item{smoothed_state_variance}{(`no_of_factors` \eqn{\times}{x} 
#'   (`no_of_factors * no_of_obs`)) matrix, where each (`no_of_factors` 
#'   \eqn{\times}{x} `no_of_factors`) block represents the smoother uncertainty
#'   at time point \eqn{t}{t}}
#'   \item{factor_var_lag_order}{Integer order of the VAR process in the state 
#'   equation.}
#'   \item{llt_success_code}{Integer indicating the status of the Cholesky
#'   factorization: `0` = LLT succeeded, `-1` = LLT failed but LDLT succeeded,
#'   `-2` = both failed and errors are treated as uncorrelated.}
#' }
#'
#' @author
#' Domenic Franjic
#' 
#' @references
#' \insertRef{koopman2000fast}{TwoStepSDFM}
#' 
#' \insertRef{eigenweb}{TwoStepSDFM}
#'
#' @examples
#' data(factor_model)
#' no_of_factors <- dim(factor_model$factors)[2]
#' factor_lag_order <- dim(factor_model$trans_var_coeff)[2] / no_of_factors
#' filter_fit <- kalmanFilterSmoother(data = factor_model$data, delay = factor_model$delay, 
#'                                    no_of_factors = no_of_factors, loading_matrix = factor_model$loading_matrix, 
#'                                    meas_error_var_cov = factor_model$meas_error_var_cov, trans_error_var_cov = factor_model$trans_error_var_cov,
#'                                    trans_var_coeff = factor_model$trans_var_coeff, factor_lag_order = factor_lag_order, 
#'                                    fcast_horizon = 5, decorr_errors = TRUE,  comp_null = 1e-15, parallel = FALSE,  jitter = 1e-8)
#' print(filter_fit)
#' filter_plots <- plot(filter_fit)
#' filter_plots$`Factor Time Series Plots`
#' 
#' @export
kalmanFilterSmoother <- function(data,
                                 delay, 
                                 no_of_factors,
                                 loading_matrix,
                                 meas_error_var_cov,
                                 trans_error_var_cov,
                                 trans_var_coeff, 
                                 factor_lag_order,
                                 fcast_horizon = 0, 
                                 decorr_errors = TRUE, 
                                 comp_null = 1e-15, 
                                 parallel = FALSE, 
                                 jitter = 1e-8
) {
  func_call <- match.call()
  
  # Misshandling of the data matrix
  if(!is.zoo(data) && !is.xts(data)){
    data_r <- try(t(as.matrix(data)), silent = TRUE)
    if (inherits(data_r, "try-error")) {
      stop(paste0("data must be a matrix, convertible to a matrix or a time-series/zoo object"))
    }
  }else{
    data_r <- try(coredata(data), silent = TRUE)
    if (inherits(data_r, "try-error")) {
      stop(paste0("data must be a matrix, convertible to a matrix or a time-series/zoo object"))
    }
  }
  if(!is.numeric(data_r)){
    stop(paste0("data has non-numeric elements."))
  }
  if(any(is.infinite(data_r))){
    stop(paste0("data cannot have (-)Inf values."))
  }
  data_r[is.na(data_r)] <- 0 # Override R NAs as they seem to not get properly parsed to C++
  
  # Mishandling of delay
  no_of_vars <- dim(data_r)[2]
  no_of_obs <- dim(data_r)[1]
  if(is.null(delay)){
    delay <- matrix(rep(0, no_of_vars), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_vars)
  }
  
  # Check for NAs in the dataset outside the ragged edges
  na_ind <- FALSE
  for(col in 1:dim(data_r)[2]){
    na_ind <- any(is.na(data_r[1:(no_of_obs - delay[col]), col]))
  }
  if(na_ind){
    warning(paste0("data has NA values. AccordingThe corresponding time points will not be considered in th estimation of the loading matrix.")) 
  }
  
  # Mishandling of trans_error_var_cov
  trans_error_var_cov <- try(as.matrix(trans_error_var_cov), silent = TRUE)
  if (inherits(trans_error_var_cov, "try-error")) {
    stop(paste0("trans_error_var_cov must be a matrix or convertible to a matrix object"))
  }
  if(dim(trans_error_var_cov)[1] != dim(trans_error_var_cov)[2]){
    stop(paste0("trans_error_var_cov must be square."))
  }
  if(dim(trans_error_var_cov)[1] != no_of_factors){
    stop(paste0("trans_error_var_cov must be of dimensions no_of_factors x no_of_factors = ", no_of_factors, "x", no_of_factors))
  }
  if(!all.equal(trans_error_var_cov, t(trans_error_var_cov), tolerance = comp_null)){
    stop(paste0("trans_error_var_cov is not symmetric at tolerance level comp_null = ", comp_null, "."))
  }
  
  # Mishandling of meas_error_var_cov
  meas_error_var_cov <- try(as.matrix(meas_error_var_cov), silent = TRUE)
  if (inherits(meas_error_var_cov, "try-error")) {
    stop(paste0("meas_error_var_cov must be a matrix or convertible to a matrix object"))
  }
  if(dim(meas_error_var_cov)[1] != dim(meas_error_var_cov)[2]){
    stop(paste0("meas_error_var_cov must be square."))
  }
  if(dim(meas_error_var_cov)[1] != no_of_vars){
    stop(paste0("meas_error_var_cov must be of dimensions no_of_vars x no_of_vars = ", no_of_vars, "x", no_of_vars))
  }
  if(!all.equal(meas_error_var_cov, t(meas_error_var_cov), tolerance = comp_null)){
    stop(paste0("meas_error_var_cov is not symmetric at tolerance level comp_null = ", comp_null, "."))
  }
  
  #M´Mishandling of loading_matrix
  loading_matrix <- try(as.matrix(loading_matrix), silent = TRUE)
  if (inherits(loading_matrix, "try-error")) {
    stop(paste0("loading_matrix must be a matrix or convertible to a matrix object"))
  }
  if(dim(loading_matrix)[1] != no_of_vars || dim(loading_matrix)[2] != no_of_factors){
    stop(paste0("loading_matrix must be of dimensions no_of_vars x no_of_factors = ", no_of_vars, "x", no_of_factors))
  }
  
  # Mishandling of trans_var_coeff
  coeff_mat_parsing_error_txt <- "trans_var_coeff must be a matrix, convertible to a matrix object, or a list of no_of_factors x no_of_factors matrices of factor_lag_order"
  if(is.list(trans_var_coeff) && !is.object(trans_var_coeff)){
    if(length(trans_var_coeff) != factor_lag_order){
      stop(paste0(coeff_mat_parsing_error_txt))
    }
    for(i in 1:length(trans_var_coeff)){
      trans_var_coeff[[i]] <- try(as.matrix(trans_var_coeff[[i]]), silent = TRUE)
      if (inherits(trans_var_coeff[[i]], "try-error")) {
        stop(paste0(coeff_mat_parsing_error_txt))
      }
      if(dim(trans_var_coeff[[i]])[1] != no_of_factors ||
         dim(trans_var_coeff[[i]])[2] != no_of_factors){
        stop(paste0(coeff_mat_parsing_error_txt))
      }
    }
  }else{
    trans_var_coeff <- try(as.matrix(trans_var_coeff), silent = TRUE)
    if (inherits(trans_var_coeff, "try-error")) {
      stop(paste0(coeff_mat_parsing_error_txt))
    }
    if(dim(trans_var_coeff)[1] != no_of_factors || dim(trans_var_coeff)[2] != no_of_factors * factor_lag_order){
      stop(paste0("loading_matrix must be of dimensions no_of_factors x (no_of_factors * factor_lag_order) = ", no_of_factors, "x", no_of_factors * factor_lag_order))
    }
  }
  
  # Mishandling of number of factors
  no_of_factors <- checkPositiveSignedInteger(no_of_factors, "no_of_factors")
  if(no_of_factors == 0){
    stop("no_of_factors must be strictly positive.")
  }
  if(no_of_factors > no_of_vars){
    stop(paste0("no_of_factors must be smaller than no_of_vars."))
  }
  
  # Mishandling of number of factor_lag_order
  factor_lag_order <- checkPositiveSignedInteger(factor_lag_order, "factor_lag_order")
  if(factor_lag_order < 0){
    stop(paste0("factor_lag_order must be non-negative."))
  }
  
  # Mishandling of number of factor_lag_order
  fcast_horizon <- checkPositiveSignedInteger(fcast_horizon, "fcast_horizon")
  if(fcast_horizon < 0){
    stop(paste0("fcast_horizon must be non-negative."))
  }
  
  # Mishandling of check_rank
  decorr_errors <- checkBoolean(decorr_errors, "decorr_errors")
  
  # Mishandling of number of comp_null
  comp_null <- checkPositiveDouble(comp_null, "comp_null")
  if(comp_null < 0){
    stop(paste0("comp_null must be non-negative."))
  }
  
  # Mishandling of check_rank
  parallel <- checkBoolean(parallel, "parallel")
  
  # Mishandling of number of comp_null
  jitter <- checkPositiveDouble(jitter, "jitter")
  if(jitter < 0){
    stop(paste0("jitter must be non-negative."))
  }
  
  result <- runUVKFS(X_in = data_r,
                     delay = delay, 
                     state_var_cov = trans_error_var_cov,
                     measurement_var_cov = meas_error_var_cov,
                     loading_matrix = loading_matrix, 
                     factor_var_coefficient_matrices = trans_var_coeff, 
                     R = no_of_factors, 
                     order = factor_lag_order,
                     fcast_horizon = fcast_horizon, 
                     decorr_errors = decorr_errors, 
                     comp_null = comp_null, 
                     parallel = parallel, 
                     jitter = jitter)
  
  # Rename the results
  names(result) <- c("filtered_state_variance", "companion_form_smoothed_factors", 
                     "smoothed_state_variance", "llt_success_code")
  result$data <- data
  result$factor_var_lag_order <- factor_lag_order
  
  # Create the non-companion-form factors and loading matrix
  result$smoothed_factors <- result$companion_form_smoothed_factors[1:no_of_factors, 1:(no_of_obs + fcast_horizon), drop = FALSE]
  rownames(result$smoothed_factors) <- paste0("Factor ", 1:no_of_factors)
  result$loading_matrix_estimate <- result$loading_matrix_estimate[, 1:no_of_factors, drop = FALSE]
  result$data <- data
  no_of_cols <- no_of_obs * no_of_factors
  block_size <- factor_lag_order * no_of_factors
  temp_smoothed_state_variance <- result$smoothed_state_variance
  result$smoothed_state_variance <- matrix(NaN, no_of_factors, no_of_factors * (no_of_obs + fcast_horizon))
  result$smoothed_state_variance[, 1:no_of_factors] <- temp_smoothed_state_variance[1:no_of_factors, 1:no_of_factors]
  for(curr_obs in 2:(no_of_obs + fcast_horizon)){
    block_starting_index <- (curr_obs - 1) * block_size + 1
    block_ending_index <- block_starting_index + no_of_factors - 1
    factor_block_starting_ind <- (curr_obs - 1) * no_of_factors + 1
    factor_block_ending_ind <- factor_block_starting_ind + no_of_factors - 1
    if(curr_obs <= no_of_obs){
      result$smoothed_state_variance[, factor_block_starting_ind:factor_block_ending_ind] <- temp_smoothed_state_variance[1:no_of_factors, block_starting_index:block_ending_index]
    }else{
      result$smoothed_state_variance[, factor_block_starting_ind:factor_block_ending_ind] <- result$filtered_state_variance[1:no_of_factors, block_starting_index:block_ending_index]
    }
  }
  
  # Re-shuffle the results objects to be in a more logical ordering
  result <- result[c("data", "smoothed_factors", "smoothed_state_variance", "factor_var_lag_order", "llt_success_code")]
  
  if(is.zoo(data) || is.xts(data)){ # Also convert factors to time series
    start_vector <- c(year(time(data)[1]), month(time(data)[1]))
    result$smoothed_factors <- as.zoo(ts(t(result$smoothed_factors), start = start_vector, frequency = 12))
    colnames(result$smoothed_factors) <- paste0("Factor ", 1:no_of_factors)
  }
  
  result$call <- func_call
  class(result) <- "KFSFit"
  return(result)
}

#' @name print.KFSFit
#' @title Generic printing function for KFSFit S3 objects
#' @description
#' Print a compact summary of a `KFSFit` object.
#'
#' @param x `KFSFit` object.
#' @param ... Additional parameters.
#'
#' @return No return value, called for side effects.
#'
#' @author
#' Domenic Franjic
#' 
#' @export
print.KFSFit <- function(x, ...) {
  simulated_time_series <- is.zoo(x$smoothed_factors)
  no_of_factors <- ifelse(simulated_time_series, dim(x$smoothed_factors)[2], dim(x$smoothed_factors)[1])
  no_of_obs <- ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2])
  cat("Simulated Dynamic Factor Model\n")
  cat("=========================================================================\n")
  cat("No. of Observations                        :", no_of_obs, "\n")
  cat("No. of Factors                             :", no_of_factors, "\n")
  cat("Factor Lag Order                           :", x$factor_var_lag_order, "\n")
  if(x$llt_success_code == -1){
    cat("Info: LLT failed. Used robust LDLT instead.\n")
  }else if(x$llt_success_code == -2){
    cat("Warning: LLT and LDLT failed. Used uncorrelated errors.\n")
  }
  cat("=========================================================================\n")
  cat("Head of the factors :\n")
  max_print <- min(5, no_of_obs)
  if(simulated_time_series){
    print(head(x$smoothed_factors, max_print))
  }else{
    print(x$smoothed_factors[, 1:max_print])
  }
  cat("Tail of the factors :\n")
  if(simulated_time_series){
    print(tail(x$smoothed_factors, max_print))
  }else{
    print(x$smoothed_factors[, (dim(x$smoothed_factors)[2] - (max_print - 1)):(dim(x$smoothed_factors)[2])])
  }
  cat("=========================================================================\n")
  invisible(x)
}

#' @name plot.KFSFit
#' @title Generic plotting function for KFSFit S3 objects
#' @description
#' Create diagnostic plots for a `KFSFit` object.
#'
#' @param x `KFSFit` object.
#' @param axis_text_size Numeric size of x- and y-axis labels. Prased to ggplot2 
#' `theme(..., text = element_text(size = axis_text_size))`.
#' @param legend_title_text_size Numeric size of x- and y-axis labels. Prased to
#' ggplot2 
#' `theme(..., legend.title = element_text(size = legend_title_text_size))`.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' A named list of `patchwork`/`ggplot` objects:
#' \describe{
#'   \item{`Factor Time Series Plots`}{`patchwork`/`ggplot` object graphing the 
#'   estimated factors over time with 95% confidence bands based on the smoother
#'   uncertainty of the Kalman Filter and Smoother.}
#' }
#'
#' @author
#' Domenic Franjic
#' 
#' @export
plot.KFSFit <- function(x,
                        axis_text_size = 20, 
                        legend_title_text_size = 20, 
                        ...) {
  out_list <- list()
  if(is.zoo(x$data)){
    series_names <- colnames(x$data)
    no_of_factors <- dim(x$smoothed_factors)[2]
    time_vector <- as.Date(time(x$smoothed_factors))
    factors <- x$smoothed_factors
  }else{
    series_names <- rownames(x$data)
    no_of_factors <- dim(x$smoothed_factors)[1]
    time_vector <- 1:dim(x$smoothed_factors)[2]
    factors <- t(x$smoothed_factors)
    factors <- as.zoo(ts(factors, start = c(1, 1), frequency = 12))
  }

  out_list$`Factor Time Series Plots` <- plotFactorEstimates(factors, x$smoothed_state_variance,
                                                             no_of_factors, axis_text_size)
  
  return(out_list)
}

