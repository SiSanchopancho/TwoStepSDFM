#' @keywords internal
#' @name importsSimFM
#' ## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @import zoo
#' @import xts
#' @import lubridate
#' @import ggplot2
#' @import stats
#' @import utils
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

#' @name simFM
#' @title Draw data from a dynamic factor model.
#' @param no_of_observations Integer number of observations.
#' @param no_of_variables Integer number of Variables.
#' @param no_of_factors Integer number of factors.
#' @param loading_matrix Numeric (no_of_variables x no_of_factors) loading matrix.
#' @param meas_error_mean Numeric vector/matrix of the means of the measurement errors.
#' @param meas_error_var_cov Numeric (no_of_factors x no_of_factors) variance-covariance matrix of the measurement errors.
#' @param trans_error_var_cov Numeric (no_of_variables x no_of_variables) variance-covariance matrix of the transition errors.
#' @param trans_var_coeff Either a list of length max_factor_lag_order with each entry a numeric (no_of_factors x no_of_factors) VAR coefficient matrix or a matrix of dimensions (no_of_factors x(no_of_factors * max_factor_lag_order)) holding the VAR coefficients of the factor VAR process in each (no_of_factors x no_of_factors) block.
#' @param factor_lag_order Integer order of the VAR process in the state equation.
#' @param delay Integer vector of delays imposed onto the end of the data in months (ragged edges)
#' @param quarterfy Logical, whether or not some of the data should be aggregated to quarterly representations.
#' @param quarterly_variable_ratio Ratio of variables ought to be quarterfied.
#' @param corr Logical, whether or not the measurement error should be randomly correlated inside the function using a random correlation matrix with off-diagonal elements governed by a beta-distribution.
#' @param beta_param Parameter of the beta-distribution governing the off-diagonal elements of the variance-covariance matrix of the measurement error.
#' @param seed 32-bit unsigned integer seed for all random processes inside the function.
#' @param burn_in Integer burn-in period of the simulated data ought to be discarded at the beginning of the data.
#' @param rescale Logical, whether or not the variance of the measurement error should be rescaled by the common component to equalise the signal-to-noise ratio.
#' @param starting_date A date type object indicating the start of the dataset. If NULL (default), the function returns matrices with observations along the second dimension (i.e., time in columns). If specified, the function treats the data as a time series, aligning it accordingly.
#' @param check_stationarity Logical, whether or not the stationarity properties of the factor VAR process should be checked.
#' @param stationarity_check_threshold Threshold of the stationarity check for when to deem an eigenvalue negative beyond the numerical error.
#' @param parallel Logical, make use of Eigen internal parallel matrix operations
#' @return Returns an `simFM` object containing the DFM parameters expressed in the factor VAR companion form.
#' The `simFM` function returns a `simFM` object which, contains the following elements:
#' \describe{
#' \item{\code{`data`}}{If `starting_date` is provided, a `zoo` matrix time series. Else, `no_of_variables` x `no_of_observations` simulated data matrix.}
#' \item{\code{`factors`}}{If `starting_date` is provided, a `zoo` matrix time series. Else, `no_of_factors` x `no_of_observations` simulated factor matrix.}
#' \item{\code{`factor_var_coeff_matrix`}}{`no_of_factors` x (`no_of_factors` * `factor_lag_order`) factor VAR coefficient matrix.}
#' \item{\code{`loading_matrix`}}{Factor loading matrix.}
#' \item{\code{`measurement_error`}}{If `starting_date` is provided, a `zoo` matrix time series. Else,`no_of_variables` x `no_of_observations` simulated measurement error matrix.}
#' \item{\code{`measurement_error_var_cov`}}{Measurement error variance covariance matrix.}
#' \item{\code{`transition_error_var_cov`}}{Factor VAR process error (transition error) variance covariance matrix.}
#' \item{\code{`frequency`}}{Vector of variable frequencies}
#' \item{\code{`delay`}}{Vector of variable delays}
#' }
#' @export
simFM <- function(no_of_observations, no_of_variables, no_of_factors, loading_matrix, 
                  meas_error_mean, meas_error_var_cov, trans_error_var_cov, trans_var_coeff, 
                  factor_lag_order, delay = NULL, quarterfy = FALSE, quarterly_variable_ratio = 0, 
                  corr = FALSE, beta_param = Inf, seed = 20022024, burn_in = 1000, 
                  rescale = TRUE, starting_date = NULL,
                  check_stationarity = FALSE, stationarity_check_threshold = 1e-5,
                  parallel = FALSE) {
  
  # Mishandling of dimensionalities
  no_of_observations <- checkPositiveSignedInteger(no_of_observations, "no_of_observations")
  if(no_of_observations == 0){
    stop("no_of_observations must be strictly positive.")
  }
  
  no_of_variables <- checkPositiveSignedInteger(no_of_variables, "no_of_variables")
  if(no_of_variables == 0){
    stop("no_of_variables must be strictly positive.")
  }
  
  no_of_factors <- checkPositiveSignedInteger(no_of_factors, "no_of_factors")
  if(no_of_factors == 0){
    stop("no_of_factors must be strictly positive.")
  }
  
  if(no_of_factors > no_of_variables){
    stop(paste0("no_of_factors must be smaller than no_of_variables."))
  }
  if(no_of_observations <= no_of_variables){
    warning(paste0("no_of_variables is bigger than no_of_observations."))
  }
  
  # Mishandling of the loading matrix
  LambdaR <- checkParameterMatrix(loading_matrix, "loading_matrix", no_of_variables, no_of_factors)
  
  # Mishandling of the measurement error mean
  muR <- checkParameterMatrix(meas_error_mean, "meas_error_mean", no_of_variables, 1)
  
  
  # Mishandling of the measurement error covariance
  SigmaR <- checkParameterMatrix(meas_error_var_cov, "meas_error_var_cov", no_of_variables, no_of_variables)
  if(!isSymmetric(SigmaR)) {
    stop(paste0("meas_error_var_cov must be symmetric."))
  }
  if(any(diag(SigmaR)< 0)){
    stop(paste0("meas_error_var_cov must not have negative values along the diagonal."))
  }
  eig <- eigen(SigmaR, only.values = TRUE)$values
  if(any(eig < -1e-8)){
    warning("meas_error_var_cov may not be positive semi-definite.")
  }
  
  # Mishandling of the transition error variance
  SR <- checkParameterMatrix(trans_error_var_cov, "trans_error_var_cov", no_of_factors, no_of_factors)
  if(!isSymmetric(SR)) {
    stop(paste0("trans_error_var_cov must be symmetric."))
  }
  if(any(diag(SR)< 0)){
    stop(paste0("trans_error_var_cov must not have negative values along the diagonal."))
  }
  eig <- eigen(SR, only.values = TRUE)$values
  if(any(eig < -1e-8)){
    warning("trans_error_var_cov may not be positive semi-definite.")
  }
  
  # Mishandling of the factor VAR lag order
  factor_lag_order <- checkPositiveSignedInteger(factor_lag_order, "factor_lag_order")
  if(factor_lag_order == 0){
    stop("factor_lag_order must be strictly positive.")
  }
  
  # Misshandling of the VAR coefficient matrix
  if (is.list(trans_var_coeff) && !is.data.frame(trans_var_coeff)) { 
    # If trans_var_coeff is provided as list: Check whether each element in 
    #   trans_var_coeff has the correct dimensions and whether there is the 
    #   correct number of matrices provided.
    
    ind <- c()
    for (i in 1:length(trans_var_coeff)) {
      s <- sum(dim(trans_var_coeff[[i]]) == c(no_of_factors, no_of_factors))
      if (s != 2) {
        ind[i] <- i
      }
    }
    if (!is.null(ind)) {
      stop("The VAR coefficient matrices in trans_var_coeff must be of dimensions (no_of_factors x no_of_factors) = (", no_of_factors, "x", no_of_factors, "). The matrices with index ", paste(ind, collapse = ", "), " are of different dimensions.")
    }
    if (length(trans_var_coeff) != factor_lag_order) {
      stop("The number of VAR coefficient matrices in trans_var_coeff must equal factor_lag_order = ", factor_lag_order, " but is ", length(trans_var_coeff), ".")
    }
    
    # Store the list elements as matrix
    PhiR <- matrix(0, no_of_factors, no_of_factors * factor_lag_order)
    for (o in 1:factor_lag_order) {
      PhiR[1:no_of_factors, ((o - 1) * no_of_factors + 1):((o - 1) * no_of_factors + no_of_factors)] <- trans_var_coeff[[o]]
    }
    PhiR <- checkParameterMatrix(PhiR, "trans_var_coeff", no_of_factors, no_of_factors * factor_lag_order)
    
  } else {
    # If trans_var_coeff is provided as matrix: Do regular dimensionality checks
    PhiR <- checkParameterMatrix(trans_var_coeff, "trans_var_coeff", no_of_factors, no_of_factors * factor_lag_order)
  }
  
  # Mishandling of dealy
  if(is.null(delay)){
    delay <- matrix(rep(0, no_of_variables), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_variables)
  }
  # Mishandling of quarterfy
  quarterfy <- checkBoolean(quarterfy, "quarterfy")
  
  # Mishandling the ratio of variables ought to be quarterfied
  quarterly_variable_ratio <- checkPositiveDouble(quarterly_variable_ratio, "quarterly_variable_ratio")
  if(quarterly_variable_ratio == 0 && quarterfy){
    warning("quarterfy is set to TRUE but quarterly_variable_ratio = 0. No variables will be aggregated to quarterly data.")
  }
  if(quarterly_variable_ratio != 0 && !quarterfy){
    warning("quarterfy is set to FALSE but quarterly_variable_ratio != 0. No variables will be aggregated to quarterly data.")
  }
  if(quarterly_variable_ratio > 1){
    warning("quarterly_variable_ratio is bigger than 1. It will be set to 1 before further use. All variables will be aggregated to quarterly data.")
    quarterly_variable_ratio <- 1
  }
  
  # Mishandling of corr
  corr <- checkBoolean(corr, "corr")
  
  # Check for mishandling of the beta parameter governing the distribution of the off-diagonal elements of the meas. var.cov.
  if(is.infinite(beta_param)){
    beta_param <- .Machine$double.xmax
  }
  beta_param <- checkPositiveDouble(beta_param, "beta_param")
  if(beta_param == .Machine$double.xmax && corr){
    warning("corr is set to TRUE but beta_param = Inf. Measurement errors will not be cross-correlated.")
  }
  if(beta_param < .Machine$double.xmax && !corr){
    warning("corr is set to FALSE but beta_param != Inf. Measurement errors will not be cross-correlated.")
  }
  if(beta_param == 0){
    warning("beta_param cannot be exactly 0. It will be jittered before further use.")
    beta_param <- 1e-15
  }
  
  # Mishandling of seed
  seed <- checkPositiveSignedInteger(seed, "seed", 33)
  
  # Mishandling of burn-in
  burn_in <- checkPositiveSignedInteger(burn_in, "burn_in")
  
  # Mishandling of rescale
  rescale <- checkBoolean(rescale, "rescale")
  
  # Misshandling of the starting date
  if(!is.null(starting_date)){
    if (length(starting_date) != 1) {
      stop("starting_date must be NULL or a single string date object or single string character object convertible to a date object.")
    }
    if(!is.Date(starting_date)){
      starting_date <- try(as.Date(starting_date), silent = TRUE)
      if (inherits(starting_date, "try-error")) {
        stop(paste0("The (pseudo) starting date for the data set must be a date object or convertible to a date object."))
      }
    }
  }
  
  # Mishandling of check_stationarity
  check_stationarity <- checkBoolean(check_stationarity, "check_stationarity")
  
  # Check for mishandling of stationarity_check_threshold
  stationarity_check_threshold <- checkPositiveDouble(stationarity_check_threshold, "stationarity_check_threshold")
  if(stationarity_check_threshold == 0){
    warning("stationarity_check_threshold should not be exactly set to zero. It will be jittered before further use.")
    stationarity_check_threshold <- 1e-15
  }
  
  # Mishandling of parallel
  parallel <- checkBoolean(parallel, "parallel")
  
  # Check whether the process will result in a stationary process
  if (check_stationarity) {
    
    if(factor_lag_order > 1){
      Comp <- matrix(0, no_of_factors * factor_lag_order, no_of_factors * factor_lag_order)
      Comp[1:no_of_factors, 1:(factor_lag_order * no_of_factors)] <- PhiR
      Comp[(no_of_factors + 1):(factor_lag_order * no_of_factors), 1:(no_of_factors * (factor_lag_order - 1))] <- diag(1, no_of_factors * (factor_lag_order - 1))
    }else{
      Comp <- PhiR
    }
    
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
  
  # Generate the data
  FM <- runStaticFM(
    T = no_of_observations, N = no_of_variables, S = SR, Lambda = LambdaR, mu_e = muR, Sigma_e = SigmaR, A = PhiR, order = factor_lag_order, quarterfy = quarterfy, corr = corr, beta_param = beta_param,
    m = quarterly_variable_ratio, seed = seed, R = no_of_factors, burn_in = burn_in, rescale = rescale, 
    parallel = parallel
  )
  
  # Clean up results
  names(FM) <- c("factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error_var_cov",
                 "transition_error_var_cov", "measurement_error", "data", "frequency")
  
  # Imposing the delay
  FM$delay <- delay
  if(any(delay != 0)){
    for( n in 1:length(delay)){
      if(delay[n] > 0){
        FM$data[n, (no_of_observations - delay[n] + 1):no_of_observations] <- NaN
        FM$measurement_error[n, (no_of_observations - delay[n] + 1):no_of_observations] <- NaN
      }
    }
  }
  
  # Reorder the results 
  FM <- FM[c("data", "factors", "factor_var_coeff_matrix", "loading_matrix", "measurement_error",
             "measurement_error_var_cov", "transition_error_var_cov", "frequency", "delay")]
  
  # Create pseudo names
  rownames(FM$data) <- paste0("Series ", 1:no_of_variables)
  rownames(FM$factors) <- paste0("Factor ", 1:no_of_factors)
  rownames(FM$measurement_error) <- paste0("Series ", 1:no_of_variables)
  
  if(is.Date(starting_date)){ # Turn all data to time series if stat_date is provided
    start_vector <-  c(year(starting_date), month(starting_date))
    FM$data <- as.zoo(ts(t(FM$data), start = start_vector, frequency = 12))
    FM$factors <- as.zoo(ts(t(FM$factors), start = start_vector, frequency = 12))
    FM$measurement_error <- as.zoo(ts(t(FM$measurement_error), start = start_vector, frequency = 12))
  }
  
  class(FM) <- "simFM"
  return(FM)
  
}

#' @name print.simFM
#' @title Generic printing function for simFM S3 objects
#' @param x `simFM` object.
#' @param ... Additional parameters for the plotting functions.
#' @export
print.simFM <- function(x, ...) {
  simulated_time_series <- is.zoo(x$data)
  no_of_factors <- ifelse(simulated_time_series, dim(x$factors)[2], dim(x$factors)[1])
  cat("Simulated Dynamic Factor Model\n")
  cat("=========================================================================\n")
  cat("No. of Observations        :", ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2]), "\n")
  cat("No. of Variables           :", ifelse(simulated_time_series, dim(x$data)[2], dim(x$data)[1]), "\n")
  cat("No. of Factors             :", no_of_factors, "\n")
  cat("Factor Lag Order           :", dim(x$factor_var_coeff_matrix)[2] / no_of_factors, "\n")
  if(any(x$frequency == 4)){
    cat("No. of Quarterly Variables :", sum(x$frequency == 4), "\n")
  }
  cat("=========================================================================\n")
  cat("Head of the factors :\n")
  if(simulated_time_series){
    print(head(x$factors, 5))
  }else{
    print(x$factors[, 1:5])
  }
  cat("Tail of the factors :\n")
  if(simulated_time_series){
    print(tail(x$factors, 5))
  }else{
    print(x$factors[, (dim(x$factors)[2] - 4):(dim(x$factors)[2])])
  }
  cat("Head of the observations :\n")
  if(simulated_time_series){
    print(head(x$data, 5))
  }else{
    print(x$data[, 1:5])
  }
  cat("Tail of the observations :\n")
  if(simulated_time_series){
    print(tail(x$data, 5))
  }else{
    print(x$data[, (dim(x$data)[2] - 4):(dim(x$data)[2])])
  }
  cat("=========================================================================\n")
  invisible(x)
}

#' @name plot.simFM
#' @title Generic plotting function for simFM S3 objects
#' @param x `simFM` object.
#' @param ... Additional parameters for the plotting functions.
#' @export
plot.simFM <- function(x, ...) {
  out_list <- list()
  if(is.zoo(x$data)){
    series_names <- colnames(x$data)
    no_of_factors <- dim(x$factors)[2]
    time_vector <- as.Date(time(x$data))
    factors <- x$factors
  }else{
    series_names <- rownames(x$data)
    no_of_factors <- dim(x$factors)[1]
    time_vector <- 1:dim(x$factors)[2]
    factors <- t(x$factors)
  }
  
  # Create plots for the factors
  pot_list_factors <- list()
  for(factor in 1:no_of_factors){
    
    current_factor <- data.frame(
      Date = time_vector,
      Value = factors[, factor],
      check.names = FALSE
    )
    
    current_line_plot <- ggplot(current_factor, aes(x = Date, y = Value)) +
      geom_line(colour = "black") +
      labs(title = paste0("Factor ", factor), y = "") +
      theme_minimal()
    
    pot_list_factors[[factor]] <- current_line_plot
  }
  out_list$`Factor Time Series Plots` <- patchwork::wrap_plots(pot_list_factors, ncol = 1)
  
  # Loading matrix heatmap plot
  lambda_df <- as.data.frame(x$loading_matrix)
  colnames(lambda_df) <- paste0("Factor ", 1:no_of_factors)
  if(dim(x$loading_matrix)[2] == 1){
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
  measurement_error_var_cov_df <- as.data.frame(x$measurement_error_var_cov)
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
  
  return(out_list)
}


