#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import zoo
#' @import xts
#' @import lubridate
#' @import ggplot2
#' @import stats
#' @import utils
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

#' @name twoStepDenseDFM
#' @title Estimate an DFM via PCA and Kalman Filter and Smoother according to Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. Journal of Econometrics, 164(1), 188-205.
#' @param data Numeric (no_of_variables x no_of_observations) matrix of data or zoo/xts object.
#' @param delay Integer vector of variable delays.
#' @param no_of_factors Integer number of factors.
#' @param max_factor_lag_order Integer max P of the VAR(P) process of the factors.
#' @param decorr_errors Logical, whether or not the errors should be decorrelated.
#' @param lag_estim_criterion Information criterion used for the estimation of the factor VAR order ("BIC", "AIC", "HIC").
#' @param comp_null Computational zero.
#' @param check_rank Logical, whether or not the rank of the variance-covariance matrix should be checked.
#' @param log Logical, whether or not output should be printed along the algorithm.
#' @param parallel Logical, whether or not to use Eigen's internal parallel matrix operations.
#' @param fcast_horizon Integer forecasting horizon for factor forecasts
#' @return 
#' The `twoStepDenseDFM` function returns an object `fit` which, contains the following elements:
#' \describe{
#' \item{\code{loading_matrix_estimate}}{A matrix of estimated loadings for each factor on each observed variable.}
#' \item{\code{smoothed_factors}}{The smoothed factor estimates.}
#' \item{\code{smoothed_state_variance}}{The estimated variance-covariance matrix of the residuals/errors.}
#' \item{\code{factor_var_lag_order}}{The order of the VAR(max_factor_lag_order) model used in the factor model process.}
#' \item{\code{error_var_cov_cholesky_factor}}{A matrix representing any transformations applied to the factors for error decorrelation, if \code{decorr_errors} is TRUE.}
#' }
#' @export
twoStepDenseDFM <- function(data,
                       delay,
                       no_of_factors,
                       max_factor_lag_order = 10,
                       decorr_errors = TRUE,
                       lag_estim_criterion = "BIC",
                       comp_null = 1e-15,
                       check_rank = FALSE,
                       log = FALSE,
                       parallel = FALSE,
                       fcast_horizon = 0) {
  
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
  
  # Mishandling of delay
  no_of_variables <- dim(data_r)[2]
  no_of_observations <- dim(data_r)[1]
  if(is.null(delay)){
    delay <- matrix(rep(0, no_of_variables), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_variables)
  }
  
  # Check for NAs in the dataset outside the ragged edges
  na_ind <- c()
  for(col in 1:dim(data_r)[2]){
    na_ind[col] <- any(is.na(data_r[1:(no_of_observations - delay[col]), col]))
  }
  if(any(na_ind)){
    stop(paste0("data has NA values outside the ragged edges for the following variables: ", which(na_ind))) 
  }
  data_r[is.na(data_r)] <- 0 # Override R NAs as they seem to not get properly parsed to C++
  
  # Misshandling of dimensions
  if(no_of_variables >= no_of_observations){
    stop(paste0("Too few observations as no_of_variables >= no_of_observations."))
  }
  
  # Mishandling of number of factors
  no_of_factors <- checkPositiveSignedInteger(no_of_factors, "no_of_factors")
  if(no_of_factors == 0){
    stop("no_of_factors must be strictly positive.")
  }
  if(no_of_factors > no_of_variables){
    stop(paste0("no_of_factors must be smaller than no_of_variables."))
  }
  
  # Mishandling of number max_factor_lag_order
  max_factor_lag_order <- checkPositiveSignedInteger(max_factor_lag_order, "max_factor_lag_order")
  if(max_factor_lag_order == 0){
    stop(paste0("max_factor_lag_order must be strictly positve."))
  }
  
  # Mishandling decorr_errors
  decorr_errors <- checkBoolean(decorr_errors, "decorr_errors")
  
  # Mishandling of lag_estim_criterion
  if(is.null(lag_estim_criterion))
  {
    stop(paste0("lag_estim_criterion must be either \"BIC\", \"AIC\", or \"HIC\"."))
  }
  if(!(lag_estim_criterion %in% c("AIC", "BIC", "HIC"))){
    stop(paste0("lag_estim_criterion must be either \"BIC\", \"AIC\", or \"HIC\"."))
  }
  
  # Mishandling of comp_null
  comp_null <- checkPositiveDouble(comp_null, "comp_null")
  if(comp_null == 0){
    warning("comp_null should not be exactly 0. It will be jittered before further use.")
    comp_null <- 1e-15
  }
  
  # Mishandling of check_rank
  check_rank <- checkBoolean(check_rank, "check_rank")
  
  # Mishandling of check_rank
  log <- checkBoolean(log, "log")
  
  # Mishandling of check_rank
  parallel <- checkBoolean(parallel, "parallel")
  
  # Mishandling of comp_null
  fcast_horizon <- checkPositiveSignedInteger(fcast_horizon, "fcast_horizon")
  
  result <- runDFMKFS(
    X_in = data_r,
    delay = delay,
    R = as.integer(no_of_factors),
    order = as.integer(max_factor_lag_order),
    decorr_errors = decorr_errors,
    crit = lag_estim_criterion,
    comp_null = comp_null,
    check_rank = check_rank,
    log = log,
    KFS_conv_crit = 1e-15, # This disables intrinsic checking of thewhether or not the filter converged, as the user will be able to make this decision in non-simulation scenarios
    parallel = parallel,
    fcast_horizon = fcast_horizon
  )
  
  # Rename the results
  names(result) <- c("loading_matrix_estimate", "filtered_state_variance", "companion_form_smoothed_factors",
                     "smoothed_state_variance", "error_var_cov_cholesky_factor",
                     "factor_var_lag_order")
  
  # Create the non-companion-form factors and loading matrix
  if(result$factor_var_lag_order == 1){
    result$smoothed_factors <- result$companion_form_smoothed_factors[, 1:(no_of_observations + fcast_horizon), drop = FALSE]
  }else{
    smoothed_factors <- matrix(NaN, no_of_factors, no_of_observations + fcast_horizon)
    for(p in (result$factor_var_lag_order - 1):1){
      smoothed_factors[, (result$factor_var_lag_order - p)] <- result$companion_form_smoothed_factors[(no_of_factors * p + 1):(no_of_factors * (p + 1)), 1, drop = FALSE]
    }
    smoothed_factors[, result$factor_var_lag_order:(no_of_observations + fcast_horizon)] <- result$companion_form_smoothed_factors[1:no_of_factors, 1:(dim(result$companion_form_smoothed_factors)[2] - 1), drop = FALSE]
    result$smoothed_factors <- smoothed_factors[, 1:(no_of_observations + fcast_horizon), drop = FALSE]
  }
  rownames(result$smoothed_factors) <- paste0("Factor ", 1:no_of_factors)
  result$loading_matrix_estimate <- result$loading_matrix_estimate[, 1:no_of_factors, drop = FALSE]
  result$data <- data
  if(fcast_horizon > 0){
    no_of_cols <- no_of_observations * no_of_factors
    block_size <- fcast_horizon * no_of_factors
    result$smoothed_state_variance <- cbind(result$smoothed_state_variance,
                                            result$filtered_state_variance[, (no_of_cols - block_size + 1):no_of_cols, drop = FALSE])
  }
  
  # Re-shuffle the results objects to be in a more logical ordering
  result <- result[c("data", "loading_matrix_estimate", "smoothed_factors", "smoothed_state_variance",
                     "factor_var_lag_order", "error_var_cov_cholesky_factor")]
  
  # Compute the cholesky of the measurement error var.-cov.
  result$error_var_cov_cholesky_factor <- tryCatch({
    solve(result$error_var_cov_cholesky_factor)
  }, error = function(e) {
    return(paste("ERROR:", conditionMessage(e)))
  })
  if (is.matrix(result$error_var_cov_cholesky_factor)) {
    result$error_var_cov_cholesky_factor[upper.tri(result$error_var_cov_cholesky_factor)] <- 0
  }
  
  if(is.zoo(data) || is.xts(data)){ # Also convert factors to time series
    start_vector <- c(year(time(data)[1]), month(time(data)[1]))
    result$smoothed_factors <- as.zoo(ts(t(result$smoothed_factors), start = start_vector, frequency = 12))
  }
  
  class(result) <- "DFMFit"
  return(result)
}

#' @method print DFMFit
#' @export
print.DFMFit <- function(x, ...) {
  simulated_time_series <- is.zoo(x$smoothed_factors)
  no_of_factors <- ifelse(simulated_time_series, dim(x$smoothed_factors)[2], dim(x$smoothed_factors)[1])
  cat("Simulated Dynamic Factor Model\n")
  cat("=========================================================================\n")
  cat("No. of Observations                        :", ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2]), "\n")
  cat("No. of Variables                           :", ifelse(simulated_time_series, dim(x$data)[2], dim(x$data)[1]), "\n")
  cat("No. of Factors                             :", no_of_factors, "\n")
  cat("Factor Lag Order                           :", x$factor_var_lag_order, "\n")
  cat("=========================================================================\n")
  cat("Head of the factors :\n")
  if(simulated_time_series){
    print(head(x$smoothed_factors, 5))
  }else{
    print(x$smoothed_factors[, 1:5])
  }
  cat("Tail of the factors :\n")
  if(simulated_time_series){
    print(tail(x$smoothed_factors, 5))
  }else{
    print(x$smoothed_factors[, (dim(x$smoothed_factors)[2] - 4):(dim(x$smoothed_factors)[2])])
  }
  cat("Head of the loading matrix :\n")
  print(head(x$loading_matrix_estimate, 5))
  cat("Tail of the loading matrix :\n")
  print(tail(x$loading_matrix_estimate, 5))
  cat("=========================================================================\n")
  invisible(x)
}

#' @method plot DFMFit
#' @export
plot.DFMFit <- function(x, ...) {
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
  
  # Create plots for the factors
  pot_list_factors <- list()
  seq_along_dates <- 1:length(as.Date(time(factors)))
  for(factor in 1:no_of_factors){
    
    correction_factor <- 1.96 * sqrt(pmax(x$smoothed_state_variance[factor, seq(1, length(as.Date(time(factors))) * no_of_factors, by = no_of_factors)], 1e-15))
    current_factor <- data.frame(
      Date = as.Date(time(factors)),
      Value = factors[seq_along_dates, factor],
      `Upper 95%-CI` = factors[seq_along_dates, factor] + correction_factor,
      `Lower 95%-CI` = factors[seq_along_dates, factor] - correction_factor,
      check.names = FALSE
    )
    
    current_line_plot <- ggplot(current_factor, aes(x = Date, y = Value)) +
      geom_ribbon(aes(ymin = `Lower 95%-CI`, ymax = `Upper 95%-CI`), fill = "#88ccee", alpha = 0.4) +
      geom_line(colour = "black") +
      labs(title = paste0("Factor ", factor), y = "") +
      theme_minimal()
    
    pot_list_factors[[factor]] <- current_line_plot
  }
  out_list$`Factor Time Series Plots` <- patchwork::wrap_plots(pot_list_factors, ncol = 1)
  
  # Loading matrix heatmap plot
  lambda_df <- as.data.frame(x$loading_matrix_estimate)
  colnames(lambda_df) <- paste0("Factor ", 1:no_of_factors)
  if(dim(x$loading_matrix_estimate)[2] == 1){
    stacked_loadings <- lambda_df
    stacked_loadings$Factor <- "Factor 1"
  }else{
    stacked_loadings <- stack(lambda_df[, ]) 
  }
  colnames(stacked_loadings) <- c("Loading", "Factor")
  stacked_loadings$Variable <- factor(rep(series_names, no_of_factors), levels = rev(series_names))
  out_list$`Loading Matrix Heatmap` <- ggplot(stacked_loadings, aes(x = Factor, y = Variable)) +
    geom_tile(data = subset(stacked_loadings, Loading != 0), aes(fill = Loading), width = 0.9, height = 0.8) +
    geom_tile(data = subset(stacked_loadings, Loading == 0), fill = "black", width = 0.9, height = 0.8) +
    scale_fill_gradient2(low = "#88ccee", high = "#117733", na.value = "#882255", mid = "#FFFFFF") +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
          strip.text.y = element_blank(),
          panel.spacing = unit(0.01, "lines"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  # Measurement error var.-cov. matrix heatmap plot
  
  if (!is.character(x$error_var_cov_cholesky_factor)) {
    measurement_error_var_cov_df <- as.data.frame(x$error_var_cov_cholesky_factor %*% t(x$error_var_cov_cholesky_factor))
    colnames(measurement_error_var_cov_df) <- series_names
    stacked_measurement_error_var_cov <- stack(measurement_error_var_cov_df[, series_names])
    colnames(stacked_measurement_error_var_cov) <- c("(Co-)Variance", "Variable")
    stacked_measurement_error_var_cov$`(Co-)Variable` <- factor(rep(series_names, length(series_names)), levels = rev(series_names))
    out_list$`Meas. Error Var.-Cov. Matrix Heatmap` <- 
      ggplot(stacked_measurement_error_var_cov, aes(x = Variable, y = `(Co-)Variable`)) +
      geom_tile(data = stacked_measurement_error_var_cov, aes(fill = `(Co-)Variance`), width = 0.8, height = 0.8) +
      scale_fill_gradient2(low = "#88ccee", high = "#117733", na.value = "#882255", mid = "#FFFFFF") +
      scale_x_discrete(expand = c(0, 0)) +
      theme_minimal() +
      theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = -90),
            strip.text.y = element_blank(),
            panel.spacing = unit(0.01, "lines"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
  }else{
    out_list$`Meas. Error Var.-Cov. Matrix Heatmap` <- x$error_var_cov_cholesky_factor
  }
  
  return(out_list)
}

