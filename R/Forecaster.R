#' @name importsHelper
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

#' @name forecastWrapper
#' @title Internal forecasting wrapper function
#' @description
#' This function is for internal use only and may change in future releases
#' without notice. Users should use `nowcast()` instead for a stable and
#' supported interface.
#' Helper function to check parameter vectors for positive signed integer values
#' @keywords internal
forecastWrapper <- function(target_variables,
                            quarterly_predictors,
                            factors,
                            target_variable_delay,
                            quarterly_delay,
                            lag_estim_criterion,
                            max_fcast_horizon,
                            max_ar_lag_order,
                            max_predictor_lag_order
)
{
  no_of_factors <- dim(factors)[1]
  no_of_target_vars <- dim(target_variables)[1]
  qtrly_predictors_missing <- is.null(quarterly_predictors)
  if(qtrly_predictors_missing){
    no_of_qrtly_vars <- 0
  }else{
    no_of_qrtly_vars <- dim(quarterly_predictors)[1]
  }
  no_of_vars  <- no_of_qrtly_vars + no_of_target_vars
  no_of_qrtly_obs <- (dim(factors)[2] - 2) / 3
  min_fcast_horizons <- ifelse(target_variable_delay == 0, 
                               1,
                               -floor(target_variable_delay / 3) + 1)
  return_object <- list()
  
  # Store all monthly predictors and the monthly factors in a single matrix
  all_qtrly_data_delay <- c(floor(c(target_variable_delay, quarterly_delay) / 3), rep(0, no_of_factors))
  
  # Start quarterfication loop according to Mariano and Murasawa #
  
  # Note: It is implicitly assumed that the stationary monthly data starts at the 
  #   second month of the first quarter. Further, it is assumed that the data set
  #   ends with an observation in the last month of the quarter. This is handled 
  #   by the higher level function calling this wrapper.
  
  all_qrtly_data <- matrix(NaN, no_of_vars + no_of_factors, no_of_qrtly_obs)
  for(t in seq(5, dim(factors)[2], 3)){
    all_qrtly_data[1:no_of_target_vars, (t - 2)/3] <- target_variables[, t]
    if(!qtrly_predictors_missing){
      all_qrtly_data[(no_of_target_vars + 1):(no_of_target_vars + no_of_qrtly_vars), 
                     (t - 2)/3] <- quarterly_predictors[, t]
    }
    all_qrtly_data[(no_of_target_vars + no_of_qrtly_vars + 1):(no_of_vars + no_of_factors), 
                   (t - 2)/3] <- rowSums(cbind(1/3 * factors[, t, drop = FALSE],
                                               2/3 * factors[, t - 1, drop = FALSE],
                                               1 * factors[, t - 2, drop = FALSE],
                                               2/3 * factors[, t - 3, drop = FALSE],
                                               1/3 * factors[, t - 4, drop = FALSE]),
                                         na.rm = TRUE)
  }
  
  # End quarterfication loop according to Mariano and Murasawa #
  
  # Start ARDL estimation loop over the target variables #
  
  fcasts <- matrix(NaN, no_of_target_vars, max_fcast_horizon - min(min_fcast_horizons) + 1)
  for(current_target in 1:no_of_target_vars){
    
    # Start ARDL estimation loop over the predictor #
    
    current_fcasts <- matrix(NaN, no_of_vars + no_of_factors, max_fcast_horizon - min_fcast_horizons[current_target] + 1)
    for(current_predictor in 1:(no_of_vars + no_of_factors)){
      if(current_target == current_predictor){
        next # Skip using the current target_variable as a single predictor as its always included as predictor
      }      
      if(all_qtrly_data_delay[current_target] < all_qtrly_data_delay[current_predictor]){
        next # Skip a predictor if it is dalyed further back compared to the target variable (we do not expect forecasting gains from using variables that are further behind then the target)
      }
      
      rel_fcast_horizons <- min_fcast_horizons[current_target]:max_fcast_horizon + all_qtrly_data_delay[current_predictor]
      for(h in rel_fcast_horizons){
        
        # Fit the model for the specific forecasting horizon
        horizon_adjustment <- which(rel_fcast_horizons == h)
        
        horizon_specific_target <- matrix(
          all_qrtly_data[current_target, 
                         (horizon_adjustment + 1):(no_of_qrtly_obs - all_qtrly_data_delay[current_target])],
          ncol = 1)
        
        horizon_specific_ar_lag <- matrix(
          all_qrtly_data[current_target, 
                         1:(no_of_qrtly_obs - all_qtrly_data_delay[current_target] - horizon_adjustment)],
          ncol = 1)
        
        horizon_specific_predictor <- matrix(
          all_qrtly_data[current_predictor, 
                         (horizon_adjustment + 1 - h):(no_of_qrtly_obs - all_qtrly_data_delay[current_target] - h)],
          ncol = 1)
        
        ardl_fit <- runARDL(horizon_specific_target,
                            horizon_specific_ar_lag,
                            horizon_specific_predictor,
                            max(max_ar_lag_order - max(h, 0), 1), 
                            max(max_predictor_lag_order - max(h, 0), 1),
                            lag_estim_criterion)
        
        # Forecast
        forecast_predictors <- matrix(1, sum(ardl_fit$optimL_lag_order) + 3, 1) # Add three for the intercept and the "contemporaenous" observations
        
        forecast_predictors[2:(ardl_fit$optimL_lag_order[1] + 2), ] <- 
          head(all_qrtly_data[current_target, 
                              (no_of_qrtly_obs - all_qtrly_data_delay[current_target]):1],
               ardl_fit$optimL_lag_order[1] + 1)
        
        forecast_predictors[(ardl_fit$optimL_lag_order[1] + 3):(ardl_fit$optimL_lag_order[1] + ardl_fit$optimL_lag_order[2] + 3), ] <- 
          head(all_qrtly_data[current_predictor, 
                              (no_of_qrtly_obs - all_qtrly_data_delay[current_predictor]):1],
               ardl_fit$optimL_lag_order[2] + 1)
        
        current_fcasts[current_predictor, which(rel_fcast_horizons == h)] <-
          matrix(ardl_fit$coefficients, nrow = 1) %*% forecast_predictors
        
      }
      
      # End loop over the forecasting horizons #
      
    }
    
    # Store the final point forecast using simple forecast averaging for each target
    if(qtrly_predictors_missing){
      rownames(current_fcasts) <- c(rownames(target_variables), rownames(factors))
    }else{
      rownames(current_fcasts) <- c(rownames(target_variables), rownames(quarterly_predictors), rownames(factors))
    }
    return_object[[current_target]] <- current_fcasts
    names(return_object)[current_target] <- paste0("Single Predictor Forecasts ", rownames(target_variables)[current_target], collapse = "")
    fcasts[current_target, (max_fcast_horizon - min(min_fcast_horizons) - length(rel_fcast_horizons) + 2):(max_fcast_horizon - min(min_fcast_horizons) + 1)] <-
      colMeans(current_fcasts, na.rm = TRUE)
    
    # Start ARDL estimation loop over the predictor #
    
  }
  
  # End ARDL estimation loop over the target variables #
  
  rownames(fcasts) <- rownames(target_variables)
  return_object[[no_of_target_vars + 1]] <- fcasts
  names(return_object)[no_of_target_vars + 1] <- "Avg. Point Forecast"
  
  return(return_object)
  
}
