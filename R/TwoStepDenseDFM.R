#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
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

#' @name twoStepDenseDFM
#' @title Two Step Dense Dynamic Factor Model Estimator.
#' @description
#' Estimate a dense dynamic factor model with measurement equation
#' \deqn{
#'   \bm{x}_t = \bm{\Lambda} \bm{f}_{t} + \bm{\xi}_t,\quad \bm{\xi}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{\xi}),
#' }
#' and transition equation
#' \deqn{
#'   \bm{f}_t = \sum_{p=0}^P\bm{\Phi}_p \bm{f}_{t-p} + \bm{\epsilon}_t,\quad \bm{\epsilon}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{f}).
#' }
#' using principal components  analysis and the Kalman Filter and 
#' Smoother according to \insertRef{Giannone2008Nowcasting}{TwoStepSDFM} and
#' \insertRef{Doz2011Two_step}{TwoStepSDFM}.
#' 
#' @param data Numeric (no_of_vars \eqn{\times}{x} no_of_obs) matrix of data or 
#' zoo/xts object sampled at the same frequency.
#' @param delay Integer vector of variable delays.
#' @param no_of_factors Integer number of factors.
#' @param max_factor_lag_order Integer maximum order of the VAR process in the 
#' transition equation.
#' @param lag_estim_criterion Information criterion used for the estimation of 
#' the factor VAR order (`"BIC"` (default), `"AIC"`, `"HIC"`).
#' @param decorr_errors Logical, whether or not the errors should be 
#' decorrelated.
#' @param comp_null Numeric computational zero.
#' @param parallel Logical, whether or not to use Eigen's internal parallel 
#' matrix operations.
#' @param fcast_horizon Integer number of additional Filter predictions into the 
#' future.
#' @param jitter Numerical jitter for stability of internal solver algorithms. 
#' The jitter is added to the diagonal entries of the variance covariance matrix 
#' of the measurement errors.
#' 
#' @details
#' The function performs a two-step estimation procedure for dense dynamic 
#' factor models as described in \insertRef{Giannone2008Nowcasting}{TwoStepSDFM}
#' and \insertRef{Doz2011Two_step}{TwoStepSDFM}. In the first step, the factor 
#' loading matrix is estimated using PCA. In the second step the latent factors
#' are estimated using the univariate representation of the Kalman Filter and
#' Smoother \insertCite{koopman2000fast}{TwoStepSDFM}.
#' 
#' With respect to the univariate representation of the Kalman filter and 
#' smoother, `decorr_errors` indicates whether the data should be decorrelated 
#' internally prior to filtering and smoothing. `jitter` is added to the 
#' diagonal elements of the measurement variance–covariance matrix. For more 
#' details, see \code{\link{kalmanFilterSmoother}}.
#' 
#' @return 
#' An object of class SDFMFit with main components:
#' \describe{
#'   \item{data}{Original data object.}
#'   \item{loading_matrix_estim}{Numeric matrix of estimated factor loadings.}
#'   \item{smoothed_factors}{Object containing the SPCA factor estimates. The 
#'   object inherits its class from data: If data is provided as `zoo`, 
#'   `factor_estim` will be a `zoo` object. If data is provided as matrix, 
#'   `factor_estim` will be a (`no_of_factors`\eqn{\times}{x}`no_of_obs` 
#'   matrix.}
#'   \item{smoothed_state_variance}{(`no_of_factors`\eqn{\times}{x}(
#'   `no_of_factors` * `no_of_obs`)) matrix, where each (`no_of_factors`
#'   \eqn{\times}{x}`no_of_factors`) block represents the smoother uncertainty 
#'   at time point\eqn{t}{t}.}
#'   \item{factor_var_lag_order}{Integer order of the VAR process in the state 
#'   equation.}
#'   \item{error_var_cov_cholesky_factor}{Numeric lower-triangular Cholesky 
#'   factor of the estimated measurement error variance–covariance matrix.}
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
#' \insertRef{Giannone2008Nowcasting}{TwoStepSDFM}
#' 
#' \insertRef{eigenweb}{TwoStepSDFM}
#' 
#' \insertRef{Doz2011Two_step}{TwoStepSDFM}
#' 
#' @seealso
#' \code{\link{sparsePCA}}: Routine for fitting estimating a sparse factor 
#' loading matrix.
#'  
#' \code{\link{kalmanFilterSmoother}}: Routine for filtering and smoothing 
#' latent factors.
#' 
#' \code{\link{twoStepSDFM}}: Two-step estimation routine for a sparse dynamic 
#' factor model.
#' 
#' @examples
#' data(factor_model)
#' no_of_vars <- dim(factor_model$data)[2]
#' no_of_factors <- dim(factor_model$factors)[2]
#' dfm_fit <- twoStepDenseDFM(data = factor_model$data, delay = factor_model$delay, no_of_factors = no_of_factors)
#' print(dfm_fit)
#' dfm_plots <- plot(dfm_fit)
#' dfm_plots$`Factor Time Series Plots`
#' dfm_plots$`Loading Matrix Heatmap`
#' dfm_plots$`Meas. Error Var.-Cov. Matrix Heatmap`
#' dfm_plots$`Meas. Error Var.-Cov. Eigenvalue Plot`
#' 
#' @export
twoStepDenseDFM <- function (data, 
                             delay, 
                             no_of_factors, 
                             max_factor_lag_order = 10, 
                             lag_estim_criterion = "BIC", 
                             decorr_errors = TRUE, 
                             comp_null = 1e-15, 
                             parallel = FALSE, 
                             fcast_horizon = 0, 
                             jitter = 1e-08) {
  func_call <- match.call()
  
  # Mishandling of data
  if (!is.zoo(data) && !is.xts(data)) {
    data_r <- try(t(as.matrix(data)), silent = TRUE)
    if (inherits(data_r, "try-error")) {
      stop(paste0("data must be a matrix, convertible to a matrix or a time-series/zoo object"))
    }
  }
  else {
    data_r <- try(coredata(data), silent = TRUE)
    if (inherits(data_r, "try-error")) {
      stop(paste0("data must be a matrix, convertible to a matrix or a time-series/zoo object"))
    }
  }
  if (!is.numeric(data_r)) {
    stop(paste0("data has non-numeric elements."))
  }
  if (any(is.infinite(data_r))) {
    stop(paste0("data cannot have (-)Inf values."))
  }
  no_of_variables <- dim(data_r)[2]
  no_of_observations <- dim(data_r)[1]
  
  # Mishandling of delay
  if (is.null(delay)) {
    delay <- matrix(rep(0, no_of_variables), ncol = 1)
  }
  else {
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_variables)
  }
  na_ind <- FALSE
  for (col in 1:dim(data_r)[2]) {
    na_ind <- any(is.na(data_r[1:(no_of_observations - delay[col]), col]))
    if (na_ind) {
      stop(paste0("data has NA values outside the ragged edges."))
    }
  }
  obs_ind <- FALSE
  for (col in 1:dim(data_r)[2]) {
    if (delay[col] > 0) {
      obs_ind <- !all(is.na(data_r[(no_of_observations - delay[col] + 1):no_of_observations, col]))
    }
    if (obs_ind) {
      stop(paste0("data has observed values inside the ragged edges."))
    }
  }
  
  # Mishandling of dimensions and other misc. parameters
  if (no_of_variables >= no_of_observations) {
    stop(paste0("Too few observations as no_of_variables >= no_of_observations."))
  }
  no_of_factors <- checkPositiveSignedInteger(no_of_factors, "no_of_factors")
  if (no_of_factors == 0) {
    stop("no_of_factors must be strictly positive.")
  }
  if (no_of_factors > no_of_variables) {
    stop(paste0("no_of_factors must be smaller than no_of_variables."))
  }
  max_factor_lag_order <- checkPositiveSignedInteger(max_factor_lag_order, "max_factor_lag_order")
  if (max_factor_lag_order == 0) {
    stop(paste0("max_factor_lag_order must be strictly positive."))
  }
  decorr_errors <- checkBoolean(decorr_errors, "decorr_errors")
  if (is.null(lag_estim_criterion)) {
    stop(paste0("lag_estim_criterion must be either \"BIC\", \"AIC\", or \"HIC\"."))
  }
  if (!(lag_estim_criterion %in% c("AIC", "BIC", "HIC"))) {
    stop(paste0("lag_estim_criterion must be either \"BIC\", \"AIC\", or \"HIC\"."))
  }
  comp_null <- checkPositiveDouble(comp_null, "comp_null")
  if (comp_null == 0) {
    warning("comp_null should not be exactly 0. It will be jittered before further use.")
    comp_null <- 1e-15
  }
  parallel <- checkBoolean(parallel, "parallel")
  fcast_horizon <- checkPositiveSignedInteger(fcast_horizon, "fcast_horizon")
  jitter <- checkPositiveDouble(jitter, "jitter")
  
  # Estimate model parameter and filter latent factors
  result <- runDFMKFS(X_in = data_r, delay = delay, R = as.integer(no_of_factors), 
                      order = as.integer(max_factor_lag_order),  decorr_errors = decorr_errors, 
                      crit = lag_estim_criterion, comp_null = comp_null, 
                      parallel = parallel, fcast_horizon = fcast_horizon, 
                      jitter = jitter)
  
  # Re-name the results
  names(result) <- c("loading_matrix_estimate", "filtered_state_variance", 
                     "companion_form_smoothed_factors", "smoothed_state_variance", 
                     "error_var_cov_cholesky_factor", "factor_var_lag_order", 
                     "llt_success_code")
  
  # Retrieve the factors and loading matrix from the companion forms
  result$smoothed_factors <- result$companion_form_smoothed_factors[1:no_of_factors, 1:(no_of_observations + fcast_horizon), drop = FALSE]
  rownames(result$smoothed_factors) <- paste0("Factor ", 1:no_of_factors)
  result$loading_matrix_estimate <- result$loading_matrix_estimate[, 1:no_of_factors, drop = FALSE]
  
  # Store the data in the return object
  result$data <- data
  
  # Retrieve the correct KFS uncertainty blocks from the companion form
  no_of_cols <- no_of_observations * no_of_factors
  block_size <- result$factor_var_lag_order * no_of_factors
  temp_smoothed_state_variance <- result$smoothed_state_variance
  result$smoothed_state_variance <- matrix(NaN, no_of_factors, no_of_factors * (no_of_observations + fcast_horizon))
  result$smoothed_state_variance[, 1:no_of_factors] <- temp_smoothed_state_variance[1:no_of_factors, 1:no_of_factors]
  for (curr_obs in 2:(no_of_observations + fcast_horizon)) {
    block_starting_index <- (curr_obs - 1) * block_size + 1
    block_ending_index <- block_starting_index + no_of_factors - 1
    factor_block_starting_ind <- (curr_obs - 1) * no_of_factors + 1
    factor_block_ending_ind <- factor_block_starting_ind + no_of_factors - 1
    if (curr_obs <= no_of_observations) {
      result$smoothed_state_variance[, factor_block_starting_ind:factor_block_ending_ind] <- 
        temp_smoothed_state_variance[1:no_of_factors, block_starting_index:block_ending_index]
    }
    else {
      result$smoothed_state_variance[, factor_block_starting_ind:factor_block_ending_ind] <- 
        result$filtered_state_variance[1:no_of_factors, block_starting_index:block_ending_index]
    }
  }
  
  # Re-shuffle the results and cut some of them for logical coherency and debloating
  result <- result[c("data", "loading_matrix_estimate", "smoothed_factors", "smoothed_state_variance", 
                     "factor_var_lag_order", "error_var_cov_cholesky_factor", "llt_success_code")]
  
  # Compute the Cholesky factor as runDFMKFS only returns the inverse of the lower triangular Cholesky factor
  result$error_var_cov_cholesky_factor <- tryCatch({
    solve(result$error_var_cov_cholesky_factor)
  }, error = function(e) {
    return(paste("ERROR:", conditionMessage(e)))
  })
  if (is.matrix(result$error_var_cov_cholesky_factor)) {
    result$error_var_cov_cholesky_factor[upper.tri(result$error_var_cov_cholesky_factor)] <- 0
  }
  
  # Turn the smoothed factors into zoo object if data is a zoo object
  if (is.zoo(data) || is.xts(data)) {
    start_vector <- c(year(time(data)[1]), month(time(data)[1]))
    result$smoothed_factors <- as.zoo(ts(t(result$smoothed_factors), start = start_vector, frequency = 12))
  }
  
  # Collect preliminary objects in the return object
  result$call <- func_call
  result$factor_fcast_horizon <- fcast_horizon
  result$data_delay <- delay
  class(result) <- "SDFMFit"
  return(result)
}