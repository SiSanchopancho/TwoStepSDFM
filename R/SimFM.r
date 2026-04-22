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

#' @name simFM
#' @title Simulate Dynamic Factor Models.
#' @description
#' Simulate data from a linear Gaussian state-space model (latent factor model),
#' with measurement equation
#' \deqn{
#'   \bm{x}_t = \bm{\Lambda} \bm{f}_{t} + \bm{\xi}_t,\quad \bm{\xi}_t \sim \mathcal{N}(\bm{\mu}, \bm{\Sigma}_{\xi}),
#' }
#' and transition equation
#' \deqn{
#'   \bm{f}_t = \sum_{p=1}^P\bm{\Phi}_p \bm{f}_{t-p} + \bm{\epsilon}_t,\quad \bm{\epsilon}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{f}).
#' }
#' for t = 1, \dots, T, as is used in, among others, 
#' \insertRef{franjic2024nowcasting}{TwoStepSDFM}.
#' 
#' @param no_of_obs Integer number of observations.
#' @param no_of_vars Integer number of Variables.
#' @param no_of_factors Integer number of factors.
#' @param loading_matrix Numeric (`no_of_vars` \eqn{\times}{x} `no_of_factors`) 
#' loading matrix.
#' @param meas_error_mean Numeric vector of the means of the measurement errors.
#' @param meas_error_var_cov Numeric (`no_of_vars` \eqn{\times}{x} 
#' `no_of_vars`) variance-covariance matrix of the measurement errors.
#' @param trans_error_var_cov Numeric (`no_of_factors` \eqn{\times}{x} 
#' `no_of_factors`) variance-covariance matrix of the transition errors.
#' @param trans_var_coeff Either a list of length `max_factor_lag_order` with 
#' each entry a numeric (`no_of_factors` \eqn{\times}{x} `no_of_factors`) VAR 
#' coefficient matrix or a matrix of dimensions (`no_of_factors` \eqn{\times}{x}
#' (`no_of_factors * max_factor_lag_order`)) holding the VAR coefficients of the
#' factor VAR process in each (`no_of_factors` \eqn{\times}{x} `no_of_factors`) 
#' block.
#' @param factor_lag_order Integer order of the VAR process in the transition 
#' equation.
#' @param delay Integer vector of delays imposed onto the end of the data 
#' (ragged edges).
#' @param quarterfy Logical, whether or not some of the data should be 
#' aggregated to quarterly representations.
#' @param quarterly_variable_ratio Ratio of variables ought to be quarterfied.
#' @param corr Logical, whether or not the measurement error should be randomly 
#' correlated inside the function using a random correlation matrix with 
#' off-diagonal elements governed by a beta-distribution.
#' @param beta_param Parameter of the beta-distribution governing the 
#' off-diagonal elements of the variance-covariance matrix of the measurement 
#' error.
#' @param seed 32-bit unsigned integer seed for all random processes inside the 
#' function.
#' @param burn_in Integer burn-in period of the simulated data ought to be 
#' discarded at the beginning of the sample.
#' @param rescale Logical, whether or not the variance of the measurement error 
#' should be rescaled by the common component to equalise the 
#' signal-to-noise ratio.
#' @param starting_date A date type object indicating the start of the dataset. 
#' If NULL (default), the function returns matrices with observations along the 
#' second dimension (i.e., time in columns). If specified, the function treats 
#' the data as a time series and returns a `zoo` object.
#' @param check_stationarity Logical, whether or not the stationarity properties
#'  of the factor VAR process should be checked.
#' @param stationarity_check_threshold Threshold of the stationarity check for 
#' when to deem an eigenvalue numerically negative.
#' @param parallel Logical, make use of Eigen internal parallel matrix 
#' operations.
#' 
#' @details
#' The `delay` vector indicates the number of observations at the end of the 
#' sample that will be set to `NA` for each variable. Here, `delay` refers to 
#' the number of months for monthly data and the number of quarters for 
#' quarterly data. For example, consider `delay <- c(1, 1)` and assume the 
#' variable with index `1` will be quarterfied. In that case, the variable with 
#' index `1` will be delayed by 1 quarter, i.e., it will be missing 3 
#' observations at the end of the panel. The variable with index `2` will be 
#' delayed by 1 month, i.e., it will be missing 1 observation at the end of the 
#' panel. This convention differs from the `delay` object of the `SimulData` 
#' class this function returns. There, `delay` represents the number of months 
#' since the most recent publication. For monthly variables, these values 
#' coincide, but for quarterly variables they are inherently different.
#' 
#' If `quarterfy = TRUE`, `floor(quarterly_variable_ratio * no_of_vars)` 
#' variables will be aggregated to a quarterly representation using the 
#' geometric mean according to 
#' \insertRef{Mariano2003new_coincident}{TwoStepSDFM}.
#' 
#' If `corr = TRUE`, the matrix `meas_error_var_cov` is internally replaced by a
#' random variance-covariance matrix: 
#' \eqn{\tilde{\bm{\Sigma}}:=\bm{S}\bm{R}\bm{S}}{new_meas_var_cov = diag_variance_sqrt * corr_matrix * diag_variance_sqrt},
#' where \eqn{\bm{S}}{diag_variance_sqrt} is a diagonal matrix with entries 
#' equal to `sqrt(diag(meas_error_var_cov))` and \eqn{\bm{R}}{corr_matrix} is a
#' random correlation matrix. \eqn{\bm{R}}{corr_matrix} is drawn according to 
#' \insertRef{lewandowski2009generating}{TwoStepSDFM} (see also
#' \url{https://stats.stackexchange.com/questions/2746/how-to-efficiently-generate-random-positive-semidefinite-correlation-matrices}).
#' The parameter `beta_param` governs the degree of cross-correlation of the 
#' off-diagonal elements. For more information see the literature cited above.
#' 
#' The random draws of the fundamental error terms are drawn within the `C++` 
#' backend. Therefore, `seed` must be provided and `set.seed()` will not 
#' guarantee reproduceability.
#' 
#' @return Returns a `SimulData` containing the following elements:
#' \describe{
#' \item{data}{If `starting_date` is provided, a `zoo` object, else, 
#' a (`no_of_vars` \eqn{\times}{x} `no_of_obs`) numeric matrix holding the 
#' simulated data.}
#' \item{factors}{If `starting_date` is provided, a `zoo` object, else a
#'  (`no_of_factors` \eqn{\times}{x} `no_of_obs`) numeric matrix holding  the 
#'  simulated latent factors.}
#' \item{trans_var_coeff}{Numeric (`no_of_factors` \eqn{\times}{x} 
#' (`no_of_factors` * `factor_lag_order`)) factor VAR coefficient matrix.}
#' \item{loading_matrix}{Numeric factor loading matrix.}
#' \item{meas_error}{If `starting_date` is provided, a `zoo` object, else 
#' a (`no_of_vars` \eqn{\times}{x} `no_of_obs`) numeric matrix holding the 
#' fundamental measurement errors.}
#' \item{meas_error_var_cov}{Numeric measurement error 
#' variance-covariance matrix.}
#' \item{trans_error_var_cov}{Numeric transition error 
#' variance-covariance matrix.}
#' \item{frequency}{Integer vector of variable frequencies.}
#' \item{delay}{Integer vector of variable delays, measured as the number of 
#' months since the latest available observation.}
#' }
#'
#' @author
#' Domenic Franjic
#' 
#' @references
#' \insertRef{Mariano2003new_coincident}{TwoStepSDFM}
#' 
#' \insertRef{lewandowski2009generating}{TwoStepSDFM}
#' 
#' \insertRef{franjic2024nowcasting}{TwoStepSDFM}
#'
#' @examples
#' seed <- 02102025
#' set.seed(seed)
#' no_of_obs <- 100
#' no_of_vars <- 50
#' no_of_factors <- 3
#' trans_error_var_cov <- diag(1, no_of_factors)
#' loading_matrix <- matrix(round(rnorm(no_of_vars * no_of_factors)), no_of_vars, no_of_factors)
#' meas_error_mean <- rep(0, no_of_vars)
#' meas_error_var_cov <- diag(1, no_of_vars)
#' trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
#' factor_lag_order <- 2
#' delay <- c(floor(rexp(no_of_vars, 1)))
#' quarterfy <- FALSE
#' quarterly_variable_ratio  <- 0
#' corr <- TRUE
#' beta_param <- 2
#' burn_in <- 999
#' starting_date <- "1970-01-01"
#' rescale <- TRUE
#' check_stationarity <- TRUE
#' stationarity_check_threshold <- 1e-10
#' factor_model <- simFM(no_of_obs = no_of_obs, no_of_vars = no_of_vars,
#'                       no_of_factors = no_of_factors, loading_matrix = loading_matrix,
#'                       meas_error_mean = meas_error_mean, 
#'                       meas_error_var_cov = meas_error_var_cov,
#'                       trans_error_var_cov = trans_error_var_cov, 
#'                       trans_var_coeff = trans_var_coeff,
#'                       factor_lag_order = factor_lag_order, delay = delay, 
#'                       quarterfy = quarterfy, 
#'                       quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr,
#'                       beta_param = beta_param, seed = seed, burn_in = burn_in, 
#'                       starting_date = starting_date, rescale = rescale, 
#'                       check_stationarity = check_stationarity,
#'                       stationarity_check_threshold = stationarity_check_threshold)
#' print(factor_model)
#' spca_plots <- plot(factor_model)
#' spca_plots$`Factor Time Series Plots`
#' spca_plots$`Loading Matrix Heatmap`
#' spca_plots$`Meas. Error Var.-Cov. Matrix Heatmap`
#' spca_plots$`Meas. Error Var.-Cov. Eigenvalue Plot`
#' spca_plots$`Data Var.-Cov. Matrix Heatmap`
#' spca_plots$`Data Var.-Cov. Eigenvalue Plot`
#'
#' @export
simFM <- function(no_of_obs, no_of_vars, no_of_factors, loading_matrix, 
                  meas_error_mean, meas_error_var_cov, trans_error_var_cov, trans_var_coeff, 
                  factor_lag_order, delay = NULL, quarterfy = FALSE, quarterly_variable_ratio = 0, 
                  corr = FALSE, beta_param = Inf, seed = 20022024, burn_in = 1000, 
                  rescale = TRUE, starting_date = NULL,
                  check_stationarity = FALSE, stationarity_check_threshold = 1e-5,
                  parallel = FALSE) {
  func_call <- match.call()
  
  # Mishandling of dimensionalities
  no_of_obs <- checkPositiveSignedInteger(no_of_obs, "no_of_obs")
  if(no_of_obs == 0){
    stop("no_of_obs must be strictly positive.")
  }
  
  no_of_vars <- checkPositiveSignedInteger(no_of_vars, "no_of_vars")
  if(no_of_vars == 0){
    stop("no_of_vars must be strictly positive.")
  }
  
  no_of_factors <- checkPositiveSignedInteger(no_of_factors, "no_of_factors")
  if(no_of_factors == 0){
    stop("no_of_factors must be strictly positive.")
  }
  
  if(no_of_factors > no_of_vars){
    stop(paste0("no_of_factors must be smaller than no_of_vars."))
  }
  if(no_of_obs <= no_of_vars){
    warning(paste0("no_of_vars is bigger than no_of_obs."))
  }
  
  # Mishandling of the loading matrix
  LambdaR <- checkParameterMatrix(loading_matrix, "loading_matrix", no_of_vars, no_of_factors)
  
  # Mishandling of the measurement error mean
  muR <- checkParameterMatrix(meas_error_mean, "meas_error_mean", no_of_vars, 1)
  
  
  # Mishandling of the measurement error covariance
  SigmaR <- checkParameterMatrix(meas_error_var_cov, "meas_error_var_cov", no_of_vars, no_of_vars)
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
    delay <- matrix(rep(0, no_of_vars), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_vars)
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
    if(quarterfy && (month(starting_date) %% 3 != 1)){
      stop(paste0("If any observations will be quarterfied, `starting_month` must refer to a date that lies in the first month of its corresponding quarter."))
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
  FM <- runStaticFM(T = no_of_obs, N = no_of_vars, S = SR, Lambda = LambdaR, mu_e = muR, 
                    Sigma_e = SigmaR, A = PhiR, order = factor_lag_order, quarterfy = quarterfy, 
                    corr = corr, beta_param = beta_param, m = quarterly_variable_ratio, 
                    seed = seed, R = no_of_factors, burn_in = burn_in, rescale = rescale,
                    parallel = parallel
                    )
  
  # Clean up results
  names(FM) <- c("factors", "trans_var_coeff", "loading_matrix", "meas_error_var_cov",
                 "trans_error_var_cov", "meas_error", "data", "frequency")
  
  # Imposing the delay
  if(any(delay != 0)){
    for( n in 1:length(delay)){
      if(FM$frequency[n] == 12){
        if(delay[n] > 0){
          FM$data[n, (no_of_obs - delay[n] + 1):no_of_obs] <- NaN
          FM$meas_error[n, (no_of_obs - delay[n] + 1):no_of_obs] <- NaN
        }
      }else if(FM$frequency[n] == 4){
        # Check which month the latest observation corresponds to
        overhang <- 0
        if(FM$data[n, no_of_obs] != FM$data[n, no_of_obs - 1] &&
           FM$data[n, no_of_obs - 1] == FM$data[n, no_of_obs - 2]){
          overhang <- 1
        }
        if(FM$data[n, no_of_obs] != FM$data[n, no_of_obs - 1] &&
           FM$data[n, no_of_obs - 1] != FM$data[n, no_of_obs - 2]){
          overhang <- 2
        }
        delay[n] <- delay[n] * 3 + overhang
        if(delay[n] > 0){
          FM$data[n, (no_of_obs - delay[n] + 1):no_of_obs] <- NaN
          FM$meas_error[n, (no_of_obs - delay[n] + 1):no_of_obs] <- NaN
        }
      }
    }
  }
  FM$delay <- delay
  
  # Reorder the results 
  FM <- FM[c("data", "factors", "trans_var_coeff", "loading_matrix", "meas_error",
             "meas_error_var_cov", "trans_error_var_cov", "frequency", "delay")]
  
  # Create pseudo names
  rownames(FM$data) <- paste0("Series ", 1:no_of_vars)
  rownames(FM$factors) <- paste0("Factor ", 1:no_of_factors)
  rownames(FM$meas_error) <- paste0("Series ", 1:no_of_vars)
  
  if(is.Date(starting_date)){ # Turn all data to time series if stat_date is provided
    start_vector <-  c(year(starting_date), month(starting_date))
    FM$data <- as.zoo(ts(t(FM$data), start = start_vector, frequency = 12))
    FM$factors <- as.zoo(ts(t(FM$factors), start = start_vector, frequency = 12))
    FM$meas_error <- as.zoo(ts(t(FM$meas_error), start = start_vector, frequency = 12))
  }
  
  FM$call <- func_call
  class(FM) <- "SimulData"
  return(FM)
  
}

#' @name print.SimulData
#' @title Generic printing function for SimulData S3 objects
#' @description
#' Print a compact summary of an `SimulData` object.
#'
#' @param x `SimulData` object.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' No return value; Prints a summary to the console.
#'
#' @author
#' Domenic Franjic
#' 
#' @export
print.SimulData <- function(x, ...) {
  simulated_time_series <- is.zoo(x$data)
  no_of_factors <- ifelse(simulated_time_series, dim(x$factors)[2], dim(x$factors)[1])
  no_of_obs <- ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2])
  cat("Simulated Dynamic Factor Model\n")
  cat("=========================================================================\n")
  cat("No. of Observations        :", ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2]), "\n")
  cat("No. of Variables           :", ifelse(simulated_time_series, dim(x$data)[2], dim(x$data)[1]), "\n")
  cat("No. of Factors             :", no_of_factors, "\n")
  cat("Factor Lag Order           :", dim(x$trans_var_coeff)[2] / no_of_factors, "\n")
  max_print <- min(5, no_of_obs)
  if(any(x$frequency == 4)){
    cat("No. of Quarterly Variables :", sum(x$frequency == 4), "\n")
  }
  cat("=========================================================================\n")
  cat("Head of the factors :\n")
  if(simulated_time_series){
    print(head(x$factors, max_print))
  }else{
    print(x$factors[, 1:max_print])
  }
  cat("Tail of the factors :\n")
  if(simulated_time_series){
    print(tail(x$factors, max_print))
  }else{
    print(x$factors[, (dim(x$factors)[2] - (max_print - 1)):(dim(x$factors)[2])])
  }
  cat("Head of the observations :\n")
  if(simulated_time_series){
    print(head(x$data, max_print))
  }else{
    print(x$data[, 1:max_print])
  }
  cat("Tail of the observations :\n")
  if(simulated_time_series){
    print(tail(x$data, max_print))
  }else{
    print(x$data[, (dim(x$data)[2] - (max_print - 1)):(dim(x$data)[2])])
  }
  cat("=========================================================================\n")
  invisible(x)
}

#' @name plot.SimulData
#' @title Generic plotting function for SimulData S3 objects
#' @description
#' Create diagnostic plots for an `SimulData` object.
#'
#' @param x `SimulData` object.
#' @param axis_text_size Numeric size of x- and y-axis labels. Prased to ggplot2
#' `theme(..., text = element_text(size = axis_text_size))`.
#' @param legend_title_text_size Numeric size of x- and y-axis labels. Prased to
#' ggplot2 
#' `theme(..., legend.title = element_text(size = legend_title_text_size))`.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' A named list of plot objects:
#' \describe{
#'   \item{`Factor Time Series Plots`}{`patchwork`/`ggplot` object showing the 
#'   simulated factors over time.}
#'   \item{`Loading Matrix Heatmap`}{`ggplot` object showing a heatmap of the 
#'   simulated factor loadings. Zeros are highlighted in black.}
#'   \item{`Meas. Error Var.-Cov. Matrix Heatmap`}{`ggplot` object showing a 
#'   heatmap of the measurement error  variance-covariance matrix.}
#'   \item{`Meas. Error Var.-Cov. Eigenvalue Plot`}{`ggplot`object showing a bar 
#'   plot of the eigenvalues of the measurement error variance-covariance 
#'   matrix.}
#'   \item{`Data Var.-Cov. Matrix Heatmap`}{`ggplot` object showing a heatmap of 
#'   the data variance-covariance matrix.}
#'   \item{`Data Var.-Cov. Eigenvalue Plot`}{`ggplot` object showing a bar plot 
#'   of the eigenvalues of the data variance-covariance matrix.}
#'   }
#' @export
plot.SimulData <- function(x, 
                           axis_text_size = 20, 
                           legend_title_text_size = 20,
                           ...) 
{
  out_list <- list()
  if (is.zoo(x$data)) {
    series_names <- colnames(x$data)
    no_of_factors <- dim(x$factors)[2]
    no_of_obs <- dim(x$factors)[1]
    time_vector <- as.Date(time(x$data))
    factors <- x$factors
    data <- x$data
  }
  else {
    series_names <- rownames(x$data)
    no_of_factors <- dim(x$factors)[1]
    time_vector <- 1:dim(x$factors)[2]
    no_of_obs <- dim(x$factors)[2]
    factors <- t(x$factors)
    data <- t(x$data)
  }
  
  mtly_data <- x$data[, which(x$frequency == 12)]
  is_ther_qtly_data <- any(x$frequency == 4)
  if (is_ther_qtly_data) {
    qtly_data <- x$data[, which(x$frequency == 4)]
  }
  mtly_data_df <- fortify.zoo(mtly_data, names = "Date")
  mtly_data_long <- stack(mtly_data_df[, -1])
  names(mtly_data_long) <- c("Value", "Series")
  mtly_data_long$Series <- factor(mtly_data_long$Series)
  mtly_data_long$Date <- rep(mtly_data_df$Date, times = dim(mtly_data_df)[2] - 1)
  mtly_data_plot <- ggplot(mtly_data_long, aes(x = Date, y = Value, colour = Series)) + 
    geom_line() + 
    labs(title = ifelse(is_ther_qtly_data, "Monthly Data", ""), y = "") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1), 
          text = element_text(size = axis_text_size), 
          legend.position = "none")
  if (is_ther_qtly_data) {
    qtly_data_df <- fortify.zoo(qtly_data, names = "Date")
    qtly_data_long <- stack(qtly_data_df[, -1])
    names(qtly_data_long) <- c("Value", "Series")
    qtly_data_long$Series <- factor(qtly_data_long$Series)
    qtly_data_long$Date <- rep(qtly_data_df$Date, times = dim(qtly_data_df)[2] - 1)
    qtly_data_plot <- ggplot(qtly_data_long, aes(x = Date, y = Value, colour = Series)) + 
      geom_line() + labs(title = "Quarterly Data", y = "") + 
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1), 
            text = element_text(size = axis_text_size), 
            legend.position = "none")
    data_plots <- list()
    data_plots$qtly_data_plot <- qtly_data_plot
    data_plots$mtly_data_plot <- mtly_data_plot
    out_list$`Data Plots` <- patchwork::wrap_plots(data_plots, ncol = 1)
  }
  else {
    out_list$`Data Plots` <- mtly_data_plot
  }
  
  out_list$`Factor Time Series Plots` <- plotFactorEstimates(factors, matrix(0, no_of_factors, no_of_factors * no_of_obs), 
                                                             no_of_factors, axis_text_size)
  
  out_list$`Loading Matrix Heatmap` <- plotLoadingHeatMap(x$loading_matrix, series_names, 
                                                          no_of_factors, axis_text_size, 
                                                          legend_title_text_size)
  
  if (is.zoo(x$data)) {
    residuals <- coredata(na.omit(x$meas_error))
  }
  else {
    residuals <- na.omit(t(x$meas_error))
  }
  measurement_error_var_cov_df <- as.data.frame(t(residuals) %*% residuals * 1/(dim(residuals)[1] - 1))
  out_list$`Meas. Error Var.-Cov. Matrix Heatmap` <- plotMeasVarCovHeatmap(measurement_error_var_cov_df, 
                                                                           series_names, 
                                                                           axis_text_size, 
                                                                           legend_title_text_size)
  
  out_list$`Meas. Error Var.-Cov. Eigenvalue Plot` <- plotMeasVarCovEigenvalues(eigen(measurement_error_var_cov_df)$values, 
                                                                                no_of_factors, axis_text_size, legend_title_text_size)
  
  if (is.zoo(x$data)) {
    data <- coredata(na.omit(x$meas_error))
  }
  else {
    data <- na.omit(t(x$meas_error))
  }
  data_var_cov_df <- as.data.frame(t(data) %*% data * 1/(dim(data)[1] - 1))
  out_list$`Data Var.-Cov. Matrix Heatmap` <- plotMeasVarCovHeatmap(data_var_cov_df, 
                                                                    series_names, 
                                                                    axis_text_size, 
                                                                    legend_title_text_size)
  
  out_list$`Data Var.-Cov. Eigenvalue Plot` <- plotMeasVarCovEigenvalues(eigen(data_var_cov_df)$values, 
                                                                         no_of_factors, axis_text_size, legend_title_text_size)
  return(out_list)
}

