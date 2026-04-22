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

#' @name nowcast
#' @title Predict Mixed-Frequency Data via Dynamic Factor Models
#' @description
#' Backcast, nowcast, and forecast quarterly target variables via a sparse/dense 
#' DFM using additional monthly data with ragged edges. Forecasts are produced 
#' using all quarterly targets and a quarterly representation of latent monthly 
#' factors \insertCite{Mariano2003new_coincident}{TwoStepSDFM}. Final 
#' predictions are computed via equally weighted forecast averaging of ARDL 
#' models \insertCite{marcellino2010factor}{TwoStepSDFM} for each of the targets 
#' and quarterfied factors.
#' 
#' @param data Numeric (no_of_vars \eqn{\times}{x} no_of_obs) matrix of data or 
#' zoo/xts object sampled at mixed frequencies (quarterly and monthly).
#' @param variables_of_interest Integer vector indicating the index of all 
#' target variables.
#' @param max_fcast_horizon Maximum forecasting horizon of all targets.
#' @param delay Integer vector of variable delays, measured as the number of 
#' months since the latest available observation.
#' @param selected Integer vector of the number of selected variables for each 
#' factor.
#' @param frequency Integer vector of frequencies of the variables in the data 
#' set (currently supported: `12` for monthly and `4` for quarterly data).
#' @param no_of_factors Integer number of factors.
#' @param sparse Logical, if `TRUE` (default) a sparse DFM is used to estimate 
#' the model parameters and latent factors (see \code{\link{twoStepSDFM}}). 
#' Else, a dense DFM is used (see \code{\link{twoStepDenseDFM}}).
#' @param max_factor_lag_order Integer maximum order of the VAR process in the 
#' transition equation.
#' @param lag_estim_criterion Information criterion used for the estimation of 
#' the factor VAR order (`"BIC"` (default), `"AIC"`, `"HIC"`).
#' @param decorr_errors Logical, whether or not the errors should be 
#' decorrelated.
#' @param ridge_penalty Numeric ridge penalty.
#' @param lasso_penalty Numeric vector, lasso penalties for each factor (set to 
#' NULL to disable as stopping criterion).
#' @param max_iterations Integer maximum number of iterations of the SPCA 
#' algorithm.
#' @param max_no_steps Integer number of LARS steps (set to NULL to disable as 
#' stopping criterion).
#' @param weights Numeric vector, weights for each variable weighing the 
#' \eqn{\ell_1}{`l_1`} size constraint.
#' @param comp_null Numeric computational zero.
#' @param spca_conv_crit Numeric conversion criterion for the SPCA algorithm.
#' @param parallel Logical, whether or not to use Eigen's internal parallel 
#' matrix operations.
#' @param max_ar_lag_order Integer maximum number of lags of the target variable
#' included in the final ARDL prediction routine.
#' @param max_predictor_lag_order Integer maximum number of lags of the 
#' predictors included in the final ARDL prediction routine.
#' @param jitter Numerical jitter for stability of internal solver algorithms. 
#' The jitter is added to the diagonal entries of the variance covariance matrix 
#' of the measurement errors.
#' @param svd_method Either `"fast"` or `"precise"`. Option `"fast"` uses 
#' Eigen's BDCSVD divide and conquer method for the computation of the singular 
#' values. Option `"precise"` (default) implements the slower, but numerically 
#' more stable JacobiSVD method. 
#' 
#' @details
#' This function serves as a prediction wrapper for the
#' \code{\link{twoStepDenseDFM}} and \code{\link{twoStepSDFM}} functions. `data`
#' should be a mixed-frequency data set. Currently, only monthly and quarterly
#' data are supported. With respect to the quarterly data, the function expects
#' the realization of the quarterly observations to occur in the last month of
#' the quarter. Indicate quarterly and monthly variables via `frequency` by
#' setting the corresponding element of `frequency` to `4` for quarterly and to
#' `12` for monthly data.
#'
#' This function is only able to compute predictions for quarterly variables.
#' To impute the ragged edges of the monthly observations, and potentially
#' compute additional predictions for the monthly variables, call `predict` on
#' the `SDFMFit` object returned by \code{\link{twoStepDenseDFM}} / 
#' \code{\link{twoStepSDFM}} (see \code{\link{predict.SDFMFit}}).
#'
#' `max_fcast_horizon` sets the maximum number of forecasts predicted starting
#' from the final observation of the data set. For each target, the number of
#' backcasts and whether or not a nowcast should be computed is determined
#' internally. This is done in such a way that every missing quarterly
#' observation of the targets is predicted.
#'
#' `max_ar_lag_order` governs the maximum number of lags of the current target
#' used to predict said target in each ARDL model. `max_predictor_lag_order`
#' governs the maximum number of lags of each additional quarterly predictor,
#' including other potential targets and the aggregated factors, used to predict
#' any given target in each ARDL model. The actual number of lags is internally
#' estimated using the BIC. Setting `max_ar_lag_order = 0` disables the use of
#' target lags in its own prediction function.
#' 
#' `sparse` toggles between a sparse DFM and a dense DFM. If `sparse = FALSE`,
#' all SPCA stopping criteria and other parameters passed to the sparse
#' estimation routine are ignored (for details on these parameters see
#' \code{\link{twoStepDenseDFM}}). Parameters governing the Kalman Filter and
#' Smoother are passed directly to \code{\link{twoStepDenseDFM}} / 
#' \code{\link{twoStepSDFM}}. For details see the corresponding help pages.
#' 
#' @return 
#' The `nowcast` function returns named list containing the following objects:
#' \describe{
#'   \item{Forecasts}{Numeric matrix of the target variables and their 
#'   respective backcasts, nowcasts, and/or forecasts.}
#'   \item{SDFM Fit}{An `SDFMFit` object holding the estimates of the model 
#'   parameters and the latent factors (see \code{\link{twoStepSDFM}} or 
#'   \code{\link{twoStepDenseDFM}}).}
#' }
#' 
#' @author
#' Domenic Franjic
#' 
#' @references
#' \insertRef{Mariano2003new_coincident}{TwoStepSDFM}
#' 
#' \insertRef{marcellino2010factor}{TwoStepSDFM}
#' 
#' \insertRef{franjic2024nowcasting}{TwoStepSDFM}
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
#' \code{\link{twoStepDenseDFM}}: Two-step estimation routine for a dense 
#' dynamic factor model.
#' 
#' @examples
#' data(mixed_freq_factor_model)
#' no_of_vars <- dim(mixed_freq_factor_model$data)[2]
#' no_of_factors <- dim(mixed_freq_factor_model$factors)[2]
#' sparse_nowcast <- nowcast(data = mixed_freq_factor_model$data, variables_of_interest = c(1, 2),
#'                           max_fcast_horizon = 4, delay = mixed_freq_factor_model$delay,
#'                           selected = rep(floor(0.5 * no_of_vars), no_of_factors), 
#'                           frequency = mixed_freq_factor_model$frequency, 
#'                           no_of_factors = no_of_factors, sparse = TRUE)
#' print(sparse_nowcast)
#' dense_nowcast <- nowcast(data = mixed_freq_factor_model$data, variables_of_interest = c(1, 2),
#'                          max_fcast_horizon = 4, delay = mixed_freq_factor_model$delay,
#'                          selected = NULL, frequency = mixed_freq_factor_model$frequency, 
#'                          no_of_factors = no_of_factors, sparse = FALSE)
#' sparse_plots <- plot(sparse_nowcast)
#' sparse_plots$`Single Pred. Fcast Density Plots Series 1`
#' 
#' @export
nowcast <- function(data,
                    variables_of_interest,
                    max_fcast_horizon,
                    delay,
                    selected,
                    frequency,
                    no_of_factors,
                    sparse = TRUE,
                    max_factor_lag_order = 10,
                    lag_estim_criterion = "BIC",
                    decorr_errors = TRUE,
                    ridge_penalty = 1e-6,
                    lasso_penalty = NULL,
                    max_iterations = 1000,
                    max_no_steps = NULL,
                    weights = NULL,
                    comp_null = 1e-15,
                    spca_conv_crit = 1e-4,
                    parallel = FALSE,
                    max_ar_lag_order = 5,
                    max_predictor_lag_order = 5,
                    jitter = 1e-8,
                    svd_method = "precise") {
  func_call <- match.call()
  
  # Misshandling
  
  # The following variables will be checked inside the twoStepSDFM function:
  #   selected, no_of_factors, max_factor_lag_order, decorr_errors,
  #   ridge_penalty, lasso_penalty, max_iterations, max_no_steps, weights, 
  #   comp_null, spca_conv_crit,
  
  # Misshandling of the data matrix
  if(!is.zoo(data) && !is.xts(data)){
    stop(paste0("data must be a time-series/zoo object"))
  }
  no_of_variables <- dim(data)[2]
  no_of_observations <- dim(data)[1]
  
  # Mishandling of variables_of_interest
  variables_of_interest <- checkPositiveSignedParameterVector(variables_of_interest, "variables_of_interest", length(variables_of_interest))
  variables_of_interest <- as.vector(variables_of_interest)
  if(length(variables_of_interest) != length(unique(variables_of_interest))){
    warning("variables_of_interest has non-unique entires. Only unique entries will be used going forward.")
    variables_of_interest <- unique(variables_of_interest)
    variables_of_interest <- checkPositiveSignedParameterVector(variables_of_interest, "variables_of_interest", length(variables_of_interest))
    variables_of_interest <- as.vector(variables_of_interest)
  }
  no_of_vois <- length(variables_of_interest)
  
  # Mishandling of delay
  if(is.null(delay)){
    delay <- matrix(rep(0, no_of_variables), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_variables)
  }
  if(any(-max_fcast_horizon >= delay[variables_of_interest])) {
    stop("`max_fcast_horizon` is too recent relative to the available missing data. Adjust delays or forecast horizon.")
  }
  
  # Check for NAs in the dataset outside the ragged edges
  na_ind <- FALSE
  for(col in 1:dim(data)[2]){
    na_ind <- any(is.na(data[1:(no_of_observations - delay[col]), col]))
    if(na_ind){
      stop(paste0("data has NA values outside the ragged edges.")) 
    }
  }
  
  
  # Check for observations in the dataset inside the ragged edges
  obs_ind <- FALSE
  for(col in 1:dim(data)[2]){
    if(delay[col] > 0){
      obs_ind <- !all(is.na(data[(no_of_observations - delay[col] + 1):no_of_observations, col]))
    }
    if(obs_ind){
      stop(paste0("data has observed values inside the ragged edges.")) 
    }
  }
  
  # Mishandling of frequency
  frequency <- checkPositiveSignedParameterVector(frequency, "frequency", no_of_variables)
  if (length(frequency) != no_of_variables || any(!(frequency %in% c(4, 12)))) {
    stop(paste0("frequency has non-conform values. Currently only values 4 (quarterly data) and 12 (monthly data) are supported."))
  }
  
  # Check for mishandling of the Variables of Interest
  if(any(frequency[variables_of_interest] != 4)){
    stop(paste0("Currently, only quarterly target variables are supported."))
  }
  
  # Mishandling max_fcast_horizon
  max_fcast_horizon <- checkPositiveSignedInteger(max_fcast_horizon, "max_fcast_horizon")
  
  # Mishandling of sparse
  sparse <- checkBoolean(sparse, "sparse")
  if(!sparse){
    warning(paste0("sparse is set to FALSE. A dense DFM is used to nowcast. All LARS-EN stopping criteria are ignored."))
  }
  
  # Mishandling max_ar_lag_order
  max_ar_lag_order <- checkPositiveSignedInteger(max_ar_lag_order , "max_ar_lag_order ")
  
  # Mishandling max_predictor_lag_order
  max_predictor_lag_order <- checkPositiveSignedInteger(max_predictor_lag_order , "max_predictor_lag_order ")
  if(max_predictor_lag_order == 0){
    stop(paste0("max_predictor_lag_order must be zero."))
  }
  
  if(month(time(data))[dim(data)[1]] %in% c(3, 6, 9, 12)){
    fcast_horizon <- 0
    target_variable_delay <- delay[variables_of_interest]
    quarterly_delay <- delay[which(frequency == 4)[-variables_of_interest]]
    effective_fcast_horizon <- max_fcast_horizon
  }else if(month(time(data))[dim(data)[1]] %in% c(1, 4, 7, 10)){
    fcast_horizon <- 2
    target_variable_delay <- delay[variables_of_interest] + 3
    quarterly_delay <- delay[which(frequency == 4)[-variables_of_interest]] + 3
    effective_fcast_horizon <- max_fcast_horizon - 1
  }else if(month(time(data))[dim(data)[1]] %in% c(2, 5, 8, 11)){
    fcast_horizon <- 1
    target_variable_delay <- delay[variables_of_interest] + 3
    quarterly_delay <- delay[which(frequency == 4)[-variables_of_interest]] + 3
    effective_fcast_horizon <- max_fcast_horizon - 1
  }
  
  if(sparse){
    SDFM_fit <- twoStepSDFM(data = data[, which(frequency == 12), drop = FALSE], delay = delay[which(frequency == 12)],
                            selected = selected, no_of_factors = no_of_factors,  
                            max_factor_lag_order = max_factor_lag_order, 
                            decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion, 
                            ridge_penalty = ridge_penalty,  lasso_penalty = lasso_penalty, 
                            max_iterations = max_iterations, max_no_steps = max_no_steps,
                            weights = weights, comp_null = comp_null, spca_conv_crit = spca_conv_crit,
                            parallel = parallel, fcast_horizon = fcast_horizon, 
                            jitter = jitter, svd_method = svd_method
    )
  }else{
    SDFM_fit <- twoStepDenseDFM(data = data[, which(frequency == 12), drop = FALSE], delay = delay[which(frequency == 12)],
                                no_of_factors = no_of_factors,  
                                max_factor_lag_order = max_factor_lag_order, 
                                decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion, 
                                comp_null = comp_null, parallel = parallel, 
                                fcast_horizon = fcast_horizon, jitter = jitter
    )
  }
  
  # Prepare the data-set for the forecasting wrapper
  factor_ts <- SDFM_fit$smoothed_factors
  column_names <- c(colnames(data)[which(frequency == 4)], paste0("Factor ", 1:no_of_factors))
  fcast_data <- merge.zoo(data[, which(frequency == 4)], factor_ts)
  colnames(fcast_data) <- column_names
  
  # If data does not start at the second month of the first quarter 
  #   available, add observations at the beginning of the panel: 
  #   This convention is necessary for the aggregation scheme inside
  #   the forecastWrapper for the quarterly data. (see Mariano, R. S., & Murasawa, 
  #   Y. (2003). A new coincident index of business cycles based on monthly and 
  #   quarterly series. Journal of applied Econometrics, 18(4), 427-443.)
  if(!(month(time(fcast_data))[1] %in% c(2, 5, 8, 11))){
    temp_times <- time(fcast_data)[1]
    temp_data <- coredata(fcast_data)
    if(month(time(fcast_data))[1] %in% c(1, 4, 7, 10)){
      temp_data <- rbind(matrix(NaN, 2, dim(temp_data)[2]), temp_data)
      temp_start <- as.Date(time(fcast_data)[1]) %m-% months(2)
      fcast_data <- as.zoo(ts(temp_data, start = c(year(temp_start), month(temp_start)), frequency = 12))
      colnames(fcast_data) <- column_names
      new_no_of_obs <- no_of_observations + 2
    }else if(month(time(fcast_data))[1] %in% c(3, 6, 9, 12)){
      temp_data <- rbind(matrix(NaN, 1, dim(temp_data)[2]), temp_data)
      temp_start <- as.Date(time(fcast_data)[1]) %m-% months(1)
      fcast_data <- as.zoo(ts(temp_data, start = c(year(temp_start), month(temp_start)), frequency = 12))
      colnames(fcast_data) <- column_names
      new_no_of_obs <- no_of_observations + 1
    }
  }
  
  # Split the data set into target variables, quarterly predictors and monthly predictors
  modified_data <- t(coredata(fcast_data))
  target_variables <- modified_data[variables_of_interest, , drop = FALSE]
  quarterly_predictor_ind <- delay[which(frequency == 4)][-variables_of_interest]
  if(length(quarterly_predictor_ind) == 0){
    quarterly_predictors <- NULL
  }else{
    quarterly_predictors <- modified_data[which(frequency == 4)[-variables_of_interest], , drop = FALSE] 
  }
  factors <- modified_data[(dim(modified_data)[1] - no_of_factors + 1):(dim(modified_data)[1]), , drop = FALSE]
  
  forecasts <- forecastWrapper(target_variables = target_variables, quarterly_predictors = quarterly_predictors,
                               factors = factors, target_variable_delay = target_variable_delay, 
                               quarterly_delay = quarterly_delay, lag_estim_criterion = lag_estim_criterion,
                               max_fcast_horizon = effective_fcast_horizon, max_ar_lag_order = max_ar_lag_order,
                               max_predictor_lag_order = max_predictor_lag_order, jitter = jitter)
  
  # Create nice result object
  result <- list()
  
  # Store the forecast results
  result[[1]] <- list()
  names(result)[1] <- "Forecasts"
  
  # Store forecasts together with the quarterly target series
  qtrly_series <- data[which(month(time(data)) %in% c(3, 6, 9, 12)), variables_of_interest, drop = FALSE]
  forecast_and_series <- matrix(NaN, dim(qtrly_series)[1] + max_fcast_horizon, 2 * length(variables_of_interest))
  colnames(forecast_and_series) <- paste0(1:(2 * length(variables_of_interest)))
  for(n in 1:length(variables_of_interest)){
    forecast_and_series[1:dim(qtrly_series)[1], 2*n - 1] <- qtrly_series[, n]
    forecast_and_series[(dim(qtrly_series)[1] - floor(delay[variables_of_interest[n]] / 3) + 1):dim(forecast_and_series)[1], 2*n] <-
      forecasts$`Avg. Point Forecast`[n, (dim(forecasts$`Avg. Point Forecast`)[2] - floor(delay[variables_of_interest[n]] / 3) - max_fcast_horizon + 1):dim(forecasts$`Avg. Point Forecast`)[2]]
    colnames(forecast_and_series)[2*n - 1] <- colnames(data)[variables_of_interest[n]]
    colnames(forecast_and_series)[2*n] <- paste0("Fcast ", colnames(data)[variables_of_interest[n]])
  }
  result$Forecasts <- as.zoo(ts(forecast_and_series, start = c(year(time(data)[1]), 
                                                               quarter(time(data)[1])), frequency = 4))
  
  result$`Single Predictor Forecasts` <- list()
  for(i in 1:length(variables_of_interest)){
    result$`Single Predictor Forecasts`[[i]] <- 
      as.zoo(ts(t(forecasts[[i]]), 
                end = c(year(time(result$Forecasts)[dim(result$Forecasts)[1]]),
                        quarter(time(result$Forecasts)[dim(result$Forecasts)[1]])),
                frequency = 4
      ))
    names(result$`Single Predictor Forecasts`)[i] <- names(forecasts)[i]
  }
  
  
  result$`SDFM Fit` <- SDFM_fit
  result$call <- func_call
  class(result) <- "SDFMnowcast"
  return(result)
}

#' @name print.SDFMnowcast
#' @title Generic print function for SDFMnowcast S3 objects
#' 
#' @param x `SDFMnowcast` object.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' No return value; Prints a summary to the console.
#' 
#' @author
#' Domenic Franjic
#' 
#' @export
print.SDFMnowcast <- function(x, ...) {
  print(tail(x$Forecasts, 10))
  invisible(x)
}


#' @name plot.SDFMnowcast
#' @title Generic plotting function for SDFMnowcast S3 objects
#' @param x `SDFMnowcast` object.
#' @param axis_text_size Numeric size of x- and y-axis labels. Prased to ggplot2 
#' `theme(..., text = element_text(size = axis_text_size))`.
#' @param ... Additional parameters for the plotting functions.
#' 
#' @return
#' A named list storing of `ggplot` objects:
#' \describe{
#'   \item{`Single Pred. Fcast Density Plots x`}{`patchwork` / `ggplot` objects
#'   graphing the distribution of forecasts generated by the predictors for
#'   each prediction (backcasts, nowcasts, forecasts) for each target,
#'   respectively. Altogether, there will be as many such objects as there are
#'   targets, with `x` replaced by the column name of the target.}
#' }
#' 
#' @author
#' Domenic Franjic
#' 
#' @export
plot.SDFMnowcast <- function(x,                         
                             axis_text_size = 20,
                             ...) {
  out_list <- list()
  
  # Single Predictor Density Plots
  absolute_fcast_date <- time(x$`SDFM Fit`$data)[dim(x$`SDFM Fit`$data)[1]]
  if(dim(x$`Single Predictor Forecasts`[[1]])[2] > 2){
    for(h in 1:(dim(x$Forecasts)[2] / 2)){
      current_single_pred_raw <- x$`Single Predictor Forecasts`[[h]]
      current_single_pred <- na.omit(t(coredata(current_single_pred_raw)))
      
      plot_list <- list()
      for(horizon in 1:dim(current_single_pred)[2]){
        data_df <- data.frame(
          Series = rownames(current_single_pred),
          Value = as.numeric(current_single_pred[, horizon])
        )
        
        data_df$Series <- factor(data_df$Series, levels = rownames(current_single_pred))
        relative_fcast_date <- as.yearqtr(time(current_single_pred_raw)[horizon])
        current_horizon <- 4 * (relative_fcast_date - as.yearqtr(absolute_fcast_date))
        max_density <- max(density(data_df$Value)$y)
        current_density_plot <- ggplot(data_df, aes(x = Value)) +
          geom_density(fill = "#88ccee", alpha = 0.6, color = "#332288") +
          geom_vline(xintercept = mean(data_df$Value), colour = "#882255",
                     lty = 1) +
          geom_text(label = "Mean", y = max_density * 0.7, 
                    x = mean(data_df$Value ) + 0.1 * sqrt(var(data_df$Value)),
                    colour = "#882255") +
          geom_vline(xintercept = median(data_df$Value), colour = "#117733",
                     lty = 2) +
          geom_text(label = "Median", y = max_density * 0.8, 
                    x = median(data_df$Value ) - 0.12 * sqrt(var(data_df$Value)),
                    colour = "#117733") + 
          geom_point(aes(y = 0), size = 2) +
          geom_text(aes(label = Series, y  = 0.00), nudge_y = max_density * 0.2,
                    size = 3.5, angle = 90, color = "black") +
          labs(title = paste0(ifelse(current_horizon >= 0,
                                     ifelse(current_horizon == 0,
                                            "Nowcast", paste0(current_horizon, "-step ahead Forecast")),
                                     paste0(-current_horizon, "-step back Backcast")),
                              " for ", relative_fcast_date),
               x = "Predicted Value", y = ""
          ) +
          theme_minimal() + 
          theme(text = element_text(size = axis_text_size))
        plot_list[[horizon]] <- current_density_plot
      }
      out_list[[h]] <- patchwork::wrap_plots(plot_list, ncol = 2)
      names(out_list)[h] <- paste0("Single Pred. Fcast Density Plots ", colnames(x$Forecasts)[2 * h - 1])
    }
  }
  
  return(out_list)
}
