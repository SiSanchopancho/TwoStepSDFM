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

#' @name imputeMonthlyData
#' @title Impute monthly data.
#' @description
#' Impute missing data in monthly datasets using matrix imputation via a DFM 
#' fit,
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
#' and \insertRef{Doz2011Two_step}{TwoStepSDFM}. For more details, see 
#' \code{\link{twoStepDenseDFM}}. The resulting estimates of the loading matrix
#' and factors are then used to predict any missing observations outside the 
#' ragged edges.
#' 
#' @return 
#' A named list containing the following objects:
#' \describe{
#'   \item{imputed_data}{Imputed dataset with class inherited from `data`.}
#'   \item{predicted_data}{In-sample model predictions with class inherited from 
#'   `data`.}
#'   \item{dfm_fit}{`SDFMFit` object of the DFM fit.}
#' }
#' 
#' @author
#' Domenic Franjic
#' 
#' @references
#' 
#' \insertRef{Giannone2008Nowcasting}{TwoStepSDFM}
#' 
#' \insertRef{Doz2011Two_step}{TwoStepSDFM}
#' 
#' @examples
#' data(factor_model)
#' no_of_vars <- dim(factor_model$data)[2]
#' no_of_factors <- dim(factor_model$factors)[2]
#' factor_model$data[1:3, ] <- NA
#' factor_model$data
#' imputed_data <- imputeMonthlyData(data = factor_model$data, delay = factor_model$delay,
#'                                   no_of_factors = no_of_factors)
#' imputed_data$imputed_data
#' 
#' @export
imputeMonthlyData <- function (data, 
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
  
  if (is.zoo(data)) {
    no_of_obs <- dim(data)[1]
    no_of_vars <- dim(data)[2]
    data <- data
  }else {
    no_of_obs <- dim(data)[2]
    no_of_vars <- dim(data)[1]
    data <- as.zoo(ts(t(data), start = c(1, 1), frequency = 12))
  }
  
  latest_na_ind <- rep(0, no_of_vars)
  for (var in 1:no_of_vars) {
    missings <- which(is.na(data[, var]))
    missing_at_start <- missings[missings < (no_of_obs -   delay[var])]
    latest_na_ind[var] <- ifelse(length(missing_at_start) > 0, 
                                 max(missing_at_start),
                                 0)
  }
  if(all(latest_na_ind == 0)){
    stop(paste0("No data to impute outside of ragged edges. For filling ragged edges call predict on an SDFMFit object."))
  }
  
  data_cut <- data[-c(1:(max(latest_na_ind) + 1)), ]
  
  dfm_fit <- twoStepDenseDFM(data_cut, delay, no_of_factors, max_factor_lag_order,
                             lag_estim_criterion, decorr_errors, comp_null,
                             parallel, fcast_horizon, jitter)
  
  if (is.zoo(data)) {
    time_vector <- as.Date(time(dfm_fit$smoothed_factors))
    factors <- dfm_fit$smoothed_factors
  }else {
    time_vector <- 1:dim(dfm_fit$smoothed_factors)[2]
    factors <- t(dfm_fit$smoothed_factors)
    factors <- as.zoo(ts(factors, start = c(1, 1), frequency = 12))
  }
  
  data_pred <- factors %*% t(dfm_fit$loading_matrix_estimate)
  for (var in 1:no_of_vars) {
    missing_ind <- which(is.na(data[, var]))
    missing_ind <- missing_ind[missing_ind < (no_of_obs -   delay[var])]
    if(length(missing_ind) != 0){
      data[missing_ind, var] <- data_pred[missing_ind, var]
    }
  }
  
  out <- list()
  if (is.zoo(data)) {
    out$imputed_data <- data
    out$predicted_data <- data_pred
  }
  else {
    out$imputed_data <- t(coredata(data))
    out$predicted_data <- t(coredata(data_pred))
  }
  
  
  out$dfm_fit <- dfm_fit
  out$call <- func_call
  return(out)
}
