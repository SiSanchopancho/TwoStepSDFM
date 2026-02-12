#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import zoo
#' @import xts
#' @import lubridate
#' @import ggplot2
#' @import stats
#' @import utils
#' @import doParallel
#' @import doSNOW
#' @import foreach
#' @import parallel
NULL

# SPDX-License-Identifier: GPL-3.0-or-later
#
#  Copyright (C) 2024 Domenic Franjic
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

#' @name crossVal
#' @title Cross-validate the hyper-parameter of the mixed-frequency SDFM
#' @param data ``Zoo``/``xts`` object.
#' @param variable_of_interest Integer value indicating the index of the single target variable
#' @param fcast_horizon Integer value indicating the target horizon
#' @param delay Integer vector of predictor delays.
#' @param frequency Integer vector of frequencies of the variables in the data set (currently supported: ``12`` for monthly and ``4`` for quarterly data)
#' @param no_of_factors Integer number of factors.
#' @param seed 32-bit positive integer for drawing the random hyper-parameter candidates
#' @param min_ridge_penalty Lower bound for the sampled ridge penalty coefficient
#' @param max_ridge_penalty Upper bound for the sampled ridge penalty coefficient
#' @param cv_repititions Integer number of predictions ought to be computed for each candidate set
#' @param cv_size Integer number of candidate sets
#' @param lasso_penalty_type Character indicating the lasso penalty type. ``"selected"`` uses the number of non-zero elements of the loading matrix. ``"penalty"`` uses the lasso size constraint directly. ``"steps"`` uses the number of steps.
#' @param min_max_penalty Vector of size tow, where the first element indicates the lower and the second element indicates the upper bound of the lasso penalty equivalent. If lasso_penalty_type is set to ``"selected"`` or ``"steps"``, both elements must be strictly positive integers.
#' @param max_factor_lag_order Integer max P of the VAR(P) process of the factors.
#' @param decorr_errors Logical, whether or not the errors should be decorrelated.
#' @param lag_estim_criterion Information criterion used for the estimation of the factor VAR order (``"BIC"``, ``"AIC"``, ``"HIC"``).
#' @param ridge_penalty Ridge penalty.
#' @param lasso_penalty Numeric vector, lasso penalties for each factor (set to ``NaN`` if not intended as stopping criterion).
#' @param max_iterations Integer maximum number of iterations.
#' @param max_no_steps Integer number of max_no_steps (used for LARS-EN as an alternative).
#' @param comp_null Computational zero.
#' @param check_rank Logical, whether or not the rank of the variance-covariance matrix should be checked.
#' @param conv_crit Conversion criterion for the SPCA algorithm.
#' @param conv_threshold Conversion criterion for the coordinate descent algorithm.
#' @param parallel Logical, whether or not to run the cross-validation loop in parallel.
#' @param no_of_cores Integer number of course to use when run in parallel
#' @param max_ar_lag_order Integer maximum number of lags of the target variable ought to be included in the nowcasting step
#' @param max_predictor_lag_order Integer maximum number of lags of the predictors ought to be included in the nowcasting step
#' @param verbose Logical, whether to print some progress tracking output to the console.
#' @return An `SDFMcrossVal` object containing the cross-validation and BIC results for all hyper-parameter combinations.
#' The `crossVal` function returns an `SDFMcrossVal` object with the following elements:
#' \describe{
#'   \item{\code{CV}}{A list with components \code{CV Results} (matrix of cross-validation errors and corresponding hyper-parameter values) and
#'     \code{Min. CV} (row of \code{CV Results} with the minimum cross-validation error).}
#'   \item{\code{BIC}}{A list with components \code{BIC Results} (matrix of BIC values and corresponding hyper-parameter values) and
#'     \code{Min. BIC} (row of \code{BIC Results} with the minimum BIC).}
#' }
#' @export
crossVal <- function(data,
                     variable_of_interest,
                     fcast_horizon,
                     delay,
                     frequency,
                     no_of_factors,
                     seed,
                     min_ridge_penalty,
                     max_ridge_penalty,
                     cv_repititions,
                     cv_size,
                     lasso_penalty_type,
                     min_max_penalty,
                     max_factor_lag_order = 10,
                     decorr_errors = TRUE,
                     lag_estim_criterion = "BIC",
                     ridge_penalty = 1e-6,
                     lasso_penalty = NULL,
                     max_iterations = 1000,
                     max_no_steps = NULL,
                     comp_null = 1e-15,
                     check_rank = FALSE,
                     conv_crit = 1e-4,
                     conv_threshold = 1e-4,
                     parallel = FALSE,
                     no_of_cores = 1,
                     max_ar_lag_order = 5,
                     max_predictor_lag_order = 5,
                     verbose = TRUE) {
  
  
  # Mishandling
  
  # Mishandling of seed
  seed <- checkPositiveSignedInteger(seed, "seed", 33)
  set.seed(seed)
  
  # Misshandling of the data matrix
  if(!is.zoo(data) && !is.xts(data)){
    stop(paste0("data must be a time-series/zoo object"))
  }
  no_of_variables <- dim(data)[2]
  no_of_mtly_variables <- sum(frequency == 12)
  no_of_observations <- dim(data)[1]
  
  # Mishandling of frequency
  frequency <- checkPositiveSignedParameterVector(frequency, "frequency", no_of_variables)
  if (length(frequency) != no_of_variables || any(!(frequency %in% c(4, 12)))) {
    stop(paste0("frequency has non-conform values. Currently only values 4 (quarterly data) and 12 (monthly data) are supported."))
  }
  
  # Mishandling of delay
  if(is.null(delay)){
    delay <- matrix(rep(0, no_of_variables), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_variables)
  }
  
  # Mishandling of variable_of_interest
  variable_of_interest <- checkPositiveSignedInteger(variable_of_interest, "variable_of_interest")
  if(variable_of_interest == 0 || variable_of_interest > no_of_variables){
    stop(paste0("variable_of_interest must be a strictly positive integer between [1, no. of variables] = [1, ", no_of_variables, "]."))
  }
  if(frequency[variable_of_interest] != 4){
    stop(paste0("variable_of_interest must correpsond to a quarterly variable. Cross-validation for monthly target series is currently not supported."))
  }
  
  # Mishandling of fcast_horizon
  fcast_horizon <- checkPositiveDouble(fcast_horizon, "fcast_horizon")
  
  # Mishandling of min_ridge_penalty
  min_ridge_penalty <- checkPositiveDouble(min_ridge_penalty, "min_ridge_penalty")
  if(min_ridge_penalty == 0){
    warning("min_ridge_penalty should not be exactly 0. It will be jittered before further use.")
    min_ridge_penalty <- 1e-15
  }
  
  # Mishandling of max_ridge_penalty
  max_ridge_penalty <- checkPositiveDouble(max_ridge_penalty, "max_ridge_penalty")
  if(max_ridge_penalty < min_ridge_penalty){
    stop(paste0("max_ridge_penalty cannot be smaller than min_ridge_penalty."))
  }
  
  # Mishandling of cv_repititions
  cv_repititions <- checkPositiveSignedInteger(cv_repititions, "cv_repititions")
  if(cv_repititions == 0){
    stop(paste0("cv_repititions must be striclty positive."))
  }
  
  # Mishandling of cv_size
  cv_size <- checkPositiveSignedInteger(cv_size, "cv_size")
  if(cv_size <= 1){
    stop(paste0("cv_size must be striclty greater 1."))
  }
  
  # Mishandling of no_of_cores
  no_of_cores <- checkPositiveSignedInteger(no_of_cores, "no_of_cores")
  if(no_of_cores == 0){
    stop(paste0("no_of_cores must be a strictly positive integer."))
  }
  if(no_of_cores > floor(parallel::detectCores() / 2)){
    warning(paste0("no_of_cores is bigger than half the number of (physical) cores, i.e., ", floor(parallel::detectCores() / 2), ". For systems with multi-thhreading, it is recommended to use at most the number of physical cores."))
  }
  if(no_of_cores > parallel::detectCores()){
    stop(paste0("no_of_cores cannot be bigger as the maxmimum number of cores, i.e., ", parallel::detectCores()))
  }
  
  # Mishandling of lasso_penalty_type
  if(!(lasso_penalty_type %in% c("steps", "selected", "penalty"))){
    stop(paste0("lasso_penalty_type must be one of \"steps\", \"selected\", or \"penalty\"."))
  }
  
  # Mishandling of min_max_penalty
  if(length(min_max_penalty) != 2){
    stop("min_max_penalty must be of length 2.")
  }
  if(lasso_penalty_type %in% "steps"){
    min_max_penalty[1] <- checkPositiveSignedInteger(min_max_penalty[1], "The first element of min_max_pealty")
    min_max_penalty[2] <- checkPositiveSignedInteger(min_max_penalty[2], "The second element of min_max_pealty")
    if(min_max_penalty[1] == 0){
      stop(paste0("If lasso_penalty_type == \"steps\", the first element cannot be zero."))
    }
  }else if(lasso_penalty_type %in% "penalty"){
    min_max_penalty[1] <- checkPositiveDouble(min_max_penalty[1], "The first element of min_max_pealty")
    min_max_penalty[2] <- checkPositiveDouble(min_max_penalty[2], "The second element of min_max_pealty")
    if(min_max_penalty[1] == 0){
      warning("The first element of min_max_penalty should not be exactly 0. It will be jittered before further use.")
      min_max_penalty[1] <- 1e-15
    }
  }else if(lasso_penalty_type %in% "selected"){
    min_max_penalty[1] <- checkPositiveSignedInteger(min_max_penalty[1], "The first element of min_max_pealty")
    min_max_penalty[2] <- checkPositiveSignedInteger(min_max_penalty[2], "The second element of min_max_pealty")
    if(min_max_penalty[1] > sum(frequency == 12)){
      warning(paste0("The first element of min_max_penalty is bigger than the number of monthly variables. It is set to the number of variables for further use."))
      min_max_penalty[1] <- sum(frequency == 12)
    }
    if(min_max_penalty[2] > sum(frequency == 12)){
      warning(paste0("The second element of min_max_penalty is bigger than the number of monthly variables. It is set to the number of variables for further use."))
      min_max_penalty[2] <- sum(frequency == 12)
    }
  }
  if(min_max_penalty[1] >= min_max_penalty[2]){
    stop(paste0("The first element of min_max_penalty must not be bigger than the second element."))
  }
  
  # Checking whether the data-sets ends with a complete quarter of observations and cropping accordingly
  if(month(time(data))[no_of_observations] %% 3 != 0){
    warning(paste0("data must end at the last month of the final quarter. data is cropped for further use."))
    if(month(time(data))[no_of_observations] %in% c(1, 4, 7, 10)){
      data <- data[1:(no_of_observations - 1), ]
      no_of_observations <- dim(data)[1]
    }
    warning(paste0("data must end at the last month of the final quarter. data is cropped for further use."))
    if(month(time(data))[no_of_observations] %in% c(2, 5, 8, 11)){
      data <- data[1:(no_of_observations - 2), ]
      no_of_observations <- dim(data)[1]
    }
  }
  
  # Randomly draw the candidates for the hyper-parameters according to which LARS-EN stopping criterion should be used
  if(lasso_penalty_type %in% "steps"){
    candidates <- matrix(NaN, cv_size, 2)
    candidates[1, ] <- c(0.0, min_max_penalty[2])
  }else if(lasso_penalty_type %in% "penalty"){
    candidates <- matrix(NaN, cv_size, 1 + no_of_factors)
    candidates[1, ] <- c(0.0, rep(0.0, no_of_factors))
  }else if(lasso_penalty_type %in% "selected"){
    candidates <- matrix(NaN, cv_size, 1 + no_of_factors)
    candidates[1, ] <- c(0.0, rep(no_of_mtly_variables, no_of_factors))
  }
  
  log_min_ridge_penalty <- log(min_ridge_penalty)
  log_max_ridge_penalty <- log(max_ridge_penalty)
  for(i in 2:cv_size){
    candidates[i, 1] <- exp(runif(1, log_min_ridge_penalty, log_max_ridge_penalty))
    if(lasso_penalty_type %in% "steps"){
      candidates[i, 2] <- floor(runif(1, min_max_penalty[1], min_max_penalty[2]))
    }else if (lasso_penalty_type %in% "penalty") {
      candidates[i, 2:(no_of_factors + 1)] <- exp(runif(no_of_factors, log(min_max_penalty[1]), log(min_max_penalty[2])))
    }else if (lasso_penalty_type %in% "selected") {
      candidates[i, 2:(no_of_factors + 1)] <- floor(runif(no_of_factors, min_max_penalty[1], min_max_penalty[2]))
    }
  }
  
  cv_results <- matrix(NaN, cv_size, 1)
  bic_results <- matrix(NaN, cv_size, 1)
  if(!parallel){
    
    # Set-up progress bar
    if (verbose){
      message("Currently validating the model hyper-parameter in series.")
      pb <- txtProgressBar(max = cv_size, style = 3)
      setTxtProgressBar(pb, 0)
    }
    
    min_cv <- .Machine$double.xmax
    min_bic <- .Machine$double.xmax
    for(h in 1:cv_size){

      current_results <- 
        nowcastSpecificationHelper(cv_repititions = cv_repititions, no_of_factors = no_of_factors, no_of_variables = no_of_variables, 
                                   no_of_observations = no_of_observations, no_of_mtly_variables = no_of_mtly_variables,
                                   lasso_penalty_type = lasso_penalty_type,
                                   data = data, variable_of_interest = variable_of_interest, 
                                   fcast_horizon = fcast_horizon, delay = delay,  
                                   candidates = candidates[h, ], frequency = frequency, 
                                   max_factor_lag_order = max_factor_lag_order, 
                                   decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                                   max_iterations = max_iterations, comp_null = comp_null, 
                                   check_rank = check_rank, conv_crit = conv_crit,
                                   conv_threshold = conv_threshold, max_ar_lag_order = max_ar_lag_order,
                                   max_predictor_lag_order = max_predictor_lag_order)
      
      bic_results[h, 1] <- current_results$bic
      
      cv_results[h, 1] <- current_results$cv
      
      if(verbose){
        setTxtProgressBar(pb, h)
      }
    }
    if(verbose){
      close(pb)
    }
    
  }else if(parallel){
    
    # Set-up progress bar
    if (verbose) {
      message("Currently validating the model hyper-parameter in parallel.")
      pb <- txtProgressBar(max = cv_size, style = 3)
      progressFunc <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progressFunc)
    } else {
      pb <- NULL
      progressFunc <- NULL
      opts <- list()
    }
    
    # Set-up parallelisation
    cl <- makeCluster(no_of_cores)
    registerDoSNOW(cl)
    global_vars_to_export <- c("nowcastSpecificationHelper", "makeRaggedEdges")
    
    h_indices <- 1:cv_size
    results <- foreach(h = h_indices, 
                       .packages = c("zoo", "xts", "TwoStepSDFM", "lubridate"), 
                       .options.snow = opts,
                       .combine = 'rbind',
                       .multicombine = TRUE,
                       .export = global_vars_to_export) %dopar% {
                         
                         current_results <- 
                           nowcastSpecificationHelper(cv_repititions = cv_repititions, no_of_factors = no_of_factors, no_of_variables = no_of_variables, 
                                                      no_of_observations = no_of_observations, no_of_mtly_variables = no_of_mtly_variables,
                                                      lasso_penalty_type = lasso_penalty_type,
                                                      data = data, variable_of_interest = variable_of_interest, 
                                                      fcast_horizon = fcast_horizon, delay = delay,  
                                                      candidates = candidates[h, ], frequency = frequency, 
                                                      max_factor_lag_order = max_factor_lag_order, 
                                                      decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                                                      max_iterations = max_iterations, comp_null = comp_null, 
                                                      check_rank = check_rank, conv_crit = conv_crit,
                                                      conv_threshold = conv_threshold, max_ar_lag_order = max_ar_lag_order,
                                                      max_predictor_lag_order = max_predictor_lag_order)
                         
                         out <- as.data.frame(matrix(NaN, 1, 2))
                         colnames(out) <- names(current_results)
                         out$cv <- current_results$cv
                         out$bic <- current_results$bic
                         out
                       }
    if(verbose && !is.null(pb)){
      close(pb)
    }
    stopCluster(cl)
    
    cv_results[, 1] <- results$cv
    bic_results[, 1] <- results$bic
    
  }
  
  cv_out <- cbind(cv_results, candidates)
  bic_out <- cbind(bic_results, candidates)
  if(lasso_penalty_type %in% "selected"){
    colnames(cv_out) <- c("CV Errors",
                          "Ridge Penalty",
                          paste0(paste0("Factor ", 1:no_of_factors), " # non-zero Loadings"))
    colnames(bic_out) <- c("BIC",
                           "Ridge Penalty",
                           paste0(paste0("Factor ", 1:no_of_factors), " # non-zero Loadings"))
    
  }else if(lasso_penalty_type %in% "penalty"){
    colnames(cv_out) <- c("CV Errors",
                          "Ridge Penalty",
                          paste0(paste0("Factor ", 1:no_of_factors), " Lasso Penalty"))
    colnames(bic_out) <- c("BIC",
                           "Ridge Penalty",
                           paste0(paste0("Factor ", 1:no_of_factors), " Lasso Penalty"))
    
  }else if(lasso_penalty_type %in% "steps"){
    colnames(cv_out) <- c("CV Errors",
                          "Ridge Penalty",
                          "Maximum No. of LARS Steps")
    colnames(bic_out) <- c("BIC",
                           "Ridge Penalty",
                           "Maximum No. of LARS Steps")
    
  }
  
  result <- list()
  result$CV <- list(
    `CV Results` = cv_out,
    `Min. CV` = cv_out[which.min(cv_out[, 1]), , drop = FALSE]
  )
  result$BIC <- list(
    `BIC Results` = bic_out,
    `Min. BIC` = bic_out[which.min(bic_out[, 1]), , drop = FALSE]
  )
  
  
  class(result) <- "SDFMcrossVal"
  return(result)
}

#' @keywords internal
nowcastSpecificationHelper <- function(cv_repititions, no_of_factors, no_of_variables,
                                       no_of_observations, no_of_mtly_variables,
                                       lasso_penalty_type, data, variable_of_interest, 
                                       fcast_horizon, delay,  candidates, frequency, 
                                       max_factor_lag_order, decorr_errors, 
                                       lag_estim_criterion, max_iterations,  comp_null, 
                                       check_rank, conv_crit, conv_threshold, 
                                       max_ar_lag_order, max_predictor_lag_order){
  
  fcast_error <- c()
  fcast_ind <- 1
  for(t in rev(seq(from = delay[variable_of_interest], by = 3, length.out = cv_repititions))){

    oos_observation <- data[no_of_observations - t, variable_of_interest]
    is_data <- makeRaggedEdges(data[1:(no_of_observations - t), , drop = FALSE], delay)
    current_no_of_obs <- dim(is_data)[1]
    if(lasso_penalty_type %in% "selected"){
      current_nowcast <- nowcast(data = is_data, variables_of_interest = variable_of_interest, 
                                 max_fcast_horizon = fcast_horizon, delay = delay, 
                                 selected = candidates[2:(no_of_factors + 1)],
                                 frequency = frequency, no_of_factors = no_of_factors, 
                                 max_factor_lag_order = max_factor_lag_order, 
                                 decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                                 ridge_penalty = candidates[1], lasso_penalty = NULL,
                                 max_iterations = max_iterations, max_no_steps = NULL,
                                 comp_null = comp_null, check_rank = check_rank, 
                                 conv_crit = conv_crit, conv_threshold = conv_threshold,
                                 log = FALSE, parallel = FALSE, max_ar_lag_order = max_ar_lag_order,
                                 max_predictor_lag_order = max_predictor_lag_order)
    }else if(lasso_penalty_type %in% "penalty"){
      current_nowcast <- nowcast(data = is_data, variables_of_interest = variable_of_interest, 
                                 max_fcast_horizon = fcast_horizon, delay = delay, 
                                 selected =   rep(no_of_mtly_variables, no_of_factors),
                                 frequency = frequency, no_of_factors = no_of_factors, 
                                 max_factor_lag_order = max_factor_lag_order, 
                                 decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                                 ridge_penalty = candidates[1], lasso_penalty = candidates[2:(no_of_factors + 1)],
                                 max_iterations = max_iterations, max_no_steps = NULL,
                                 comp_null = comp_null, check_rank = check_rank, 
                                 conv_crit = conv_crit, conv_threshold = conv_threshold,
                                 log = FALSE, parallel = FALSE, max_ar_lag_order = max_ar_lag_order,
                                 max_predictor_lag_order = max_predictor_lag_order)
    }else if(lasso_penalty_type %in% "steps"){
      current_nowcast <- nowcast(data = is_data, variables_of_interest = variable_of_interest, 
                                 max_fcast_horizon = fcast_horizon, delay = delay, 
                                 selected =   rep(no_of_mtly_variables, no_of_factors),
                                 frequency = frequency, no_of_factors = no_of_factors, 
                                 max_factor_lag_order = max_factor_lag_order, 
                                 decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
                                 ridge_penalty = candidates[1], lasso_penalty = NULL,
                                 max_iterations = max_iterations, 
                                 max_no_steps = candidates[2],
                                 comp_null = comp_null, check_rank = check_rank, 
                                 conv_crit = conv_crit, conv_threshold = conv_threshold,
                                 log = FALSE, parallel = FALSE, max_ar_lag_order = max_ar_lag_order,
                                 max_predictor_lag_order = max_predictor_lag_order)
      
    }

    nowcast_indicator <- which(as.yearqtr(time(current_nowcast$Forecasts)) == as.yearqtr(time(is_data)[current_no_of_obs]))
    fcast_error[fcast_ind] <- coredata(current_nowcast$Forecasts[fcast_horizon + nowcast_indicator, 2]) - coredata(oos_observation)
    fcast_ind <- fcast_ind + 1
  }
  
  bic_h <- (mean((t(coredata(is_data[, which(frequency == 12)]))
                  - current_nowcast$`SDFM Fit`$loading_matrix_estimate
                  %*% t(coredata(current_nowcast$`SDFM Fit`$smoothed_factors[1:current_no_of_obs, , drop = FALSE]))
  )^2, na.rm = TRUE)
  + sum(current_nowcast$`SDFM Fit`$loading_matrix_estimate != 0)
  * log(no_of_variables * current_no_of_obs) / (no_of_variables * current_no_of_obs)
  )
  
  cv_h <- mean(fcast_error^2)
  
  return(list(cv = cv_h, bic = bic_h))
}

#' @method print SDFMcrossVal
#' @return No return value, called for side effects.
#' @export
print.SDFMcrossVal <- function(x, ...) {
  
  # Extrcat which LARS stopping criterion has been used
  if(any(grepl("Lasso", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    crit_label <- "Optimum Lasso Penalties"
  }else if(any(grepl("# non-zero Loadings", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    crit_label <- "# of non-zero Loadings per Factor"
  }
  
  cat("Cross-Validation Results\n")
  cat("=========================================================================\n")
  cat("Cross-Validation Error: ", x$CV$`Min. CV`[1], "\n")
  cat("Optimum Ridge Penalty : ", x$CV$`Min. CV`[2], "\n")
  if(any(grepl("Maximum No. of LARS Steps", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    cat("# of LARS steps       : ", x$CV$`Min. CV`[3], "\n")
  }else{
    cat("\n")
    cat(crit_label, "\n")
    for(n in 3:length(x$CV$`Min. CV`)){
      cat("Factor ", n - 2, "            :", x$CV$`Min. CV`[n], "\n")  
    }
  }
  cat("=========================================================================\n")
  cat("\n")
  cat("BIC Results\n")
  cat("=========================================================================\n")
  cat("BIC                   : ", x$BIC$`Min. BIC`[1], "\n")
  cat("Optimum Ridge Penalty : ", x$BIC$`Min. BIC`[2], "\n")
  if(any(grepl("Maximum No. of LARS Steps", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    cat("# of LARS steps       : ", x$BIC$`Min. BIC`[3], "\n")
  }else{
    cat("\n")
    cat(crit_label, "\n")
    for(n in 3:length(x$BIC$`Min. BIC`)){
      cat("Factor ", n - 2, "            :", x$BIC$`Min. BIC`[n], "\n")  
    }
  }
  cat("=========================================================================\n")
}

#' @method plot SDFMcrossVal
#' @return A named list of two ggplot objects: \code{CV Results} and \code{BIC Results}.
#' @export
plot.SDFMcrossVal <- function(x, ...) {
  out_list <- list()
  
  # Plot depending on which stopping criterion has been used
  if(any(grepl("Lasso", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    
    # Cross-validation results
    cv_data <- data.frame(x$CV$`CV Results`, check.names = FALSE)
    cv_data$`Lasso Penalties` <- c(NaN)
    for(i in 1:dim(cv_data)[1]){
      cv_data$`Lasso Penalties`[i] <- paste0("(", 
                                             paste0(sprintf("%.2f", cv_data[i, 3:(dim(cv_data)[2] - 1)]), collapse = ";"), 
                                             ")")
    }
    avg_lasso_penalty <- rowMeans(cv_data[, 3:(dim(cv_data)[2] - 1), drop = FALSE])
    breaks <- seq(from = min(avg_lasso_penalty, na.rm = TRUE), 
                  to = max(avg_lasso_penalty, na.rm = TRUE), 
                  length.out = 6)
    breaks[length(breaks)] <- breaks[length(breaks)] + 0.000006
    labels <- paste0("[", sprintf("%.2f", floor(breaks[1:5] / 0.01) * 0.01),
                     ", ",
                     sprintf("%.2f", floor(breaks[2:6] / 0.01) * 0.01),
                     ")")
    labels[5] <- paste0("[", sprintf("%.2f", floor(breaks[5] / 0.01) * 0.01),
                        ", ",
                        sprintf("%.2f", floor(breaks[6] / 0.01) * 0.01),
                        "]")
    binned_data_equal_width <- cut(avg_lasso_penalty,
                                   breaks = breaks,
                                   right = FALSE,
                                   include.lowest = TRUE,
                                   labels = labels)
    cv_data$`Avg. Lasso Penalty` <- as.factor(binned_data_equal_width)
    best_combo <- cv_data$`Lasso Penalties`[which.min(cv_data$`CV Errors`)]
    best_ridge <- cv_data$`Ridge Penalty`[which.min(cv_data$`CV Errors`)]
    best_cv_error <- min(cv_data$`CV Errors`)
    y_min_limit <- best_cv_error
    y_max_limit <- max(cv_data$`CV Errors`)
    out_list$`CV Results` <- ggplot(cv_data, aes(x = `Ridge Penalty`, y = `CV Errors`, colour = `Avg. Lasso Penalty`, 
                                                 shape = `Avg. Lasso Penalty`)) +
      geom_point(size = 2) +
      geom_hline(yintercept = cv_data$`CV Errors`[1], colour = "black") +
      scale_colour_manual(values =  c("#88CCEE", "#44799E", "#000000", "#41784A", "#117733"),
                          name = "Avg. Lasso Penalty") +
      scale_shape_discrete(name = "Avg. Lasso Penalty") +
      geom_point(data = subset(cv_data, `CV Errors` == min(`CV Errors`)), aes(x = `Ridge Penalty`, y = `CV Errors`), 
                 colour = "black", fill = "#882255", size = 3, shape = 22) +
      scale_y_continuous(trans = "log10", limits = c(y_min_limit, y_max_limit)) + 
      scale_x_continuous(trans = "log10") +
      annotate("text",  x = best_ridge, y = best_cv_error,
               label = best_combo, angle = 0, vjust = 1.6, hjust = 1, size = 4, color = "darkred") +
      labs(x = "log Ridge Penalty",
           y = "log CV Error") +
      theme_minimal()
    
    # BIC results
    bic_data <- data.frame(x$BIC$`BIC Results`, check.names = FALSE)
    bic_data$`Lasso Penalties` <- c(NaN)
    for(i in 1:dim(bic_data)[1]){
      bic_data$`Lasso Penalties`[i] <- paste0("(", 
                                              paste0(sprintf("%.2f", bic_data[i, 3:(dim(bic_data)[2] - 1)]), collapse = ";"), 
                                              ")")
    }
    bic_data$`Avg. Lasso Penalty` <- as.factor(binned_data_equal_width)
    best_bic_combo <- bic_data$`Lasso Penalties`[which.min(bic_data$`BIC`)]
    best_bic_ridge <- bic_data$`Ridge Penalty`[which.min(bic_data$`BIC`)]
    best_bic <- min(bic_data$`BIC`)
    y_bic_min_limit <- best_cv_error
    y_bic_max_limit <- max(bic_data$`BIC`)
    out_list$`BIC Results` <- ggplot(bic_data, aes(x = `Ridge Penalty`, y = `BIC`, colour = `Avg. Lasso Penalty`, 
                                                   shape = `Avg. Lasso Penalty`)) +
      geom_point(size = 2) +
      geom_hline(yintercept = bic_data$BIC[1], colour = "black") +
      scale_colour_manual(values =  c("#88CCEE", "#44799E", "#000000", "#41784A", "#117733"),
                          name = "Avg. Lasso Penalty") +
      scale_shape_discrete(name = "Avg. Lasso Penalty") +
      geom_point(data = subset(bic_data, `BIC` == min(`BIC`)), aes(x = `Ridge Penalty`, y = `BIC`), 
                 colour = "black", fill = "#882255", size = 3, shape = 22) +
      scale_y_continuous(limits = c(y_bic_min_limit, y_bic_max_limit)) + 
      scale_x_continuous(trans = "log10") +
      annotate("text",  x = best_bic_ridge, y = best_bic,
               label = best_bic_combo, angle = 0, vjust = 1.6, hjust = 1, size = 4, color = "darkred") +
      labs(x = "log Ridge Penalty",
           y = "BIC") +
      theme_minimal()
    
  }else if(any(grepl("# non-zero Loadings", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    
    cv_data <- data.frame(x$CV$`CV Results`, check.names = FALSE)
    cv_data$`# non-zero Loadings` <- c(NaN)
    for(i in 1:dim(cv_data)[1]){
      cv_data$`# non-zero Loadings`[i] <- paste0("(", 
                                                 paste0(cv_data[i, 3:(dim(cv_data)[2] - 1)], collapse = ";"), 
                                                 ")")
    }
    sparsity_ratios <- 1 - rowSums(cv_data[, 3:(dim(cv_data)[2] - 1), drop = FALSE]) / rowSums(cv_data[1, 3:(dim(cv_data)[2] - 1), drop = FALSE])
    breaks <- seq(from = min(sparsity_ratios, na.rm = TRUE), 
                  to = max(sparsity_ratios, na.rm = TRUE), 
                  length.out = 6)
    breaks[length(breaks)] <- breaks[length(breaks)] + 0.000006
    labels <- paste0("[", sprintf("%.2f", floor(breaks[1:5] / 0.01) * 0.01),
                     ", ",
                     sprintf("%.2f", floor(breaks[2:6] / 0.01) * 0.01),
                     ")")
    labels[5] <- paste0("[", sprintf("%.2f", floor(breaks[5] / 0.01) * 0.01),
                        ", ",
                        sprintf("%.2f", floor(breaks[6] / 0.01) * 0.01),
                        "]")
    binned_data_equal_width <- cut(sparsity_ratios,
                                   breaks = breaks,
                                   right = FALSE,
                                   include.lowest = TRUE,
                                   labels = labels)
    cv_data$`Sparsity Ratio` <- as.factor(binned_data_equal_width)
    best_combo <- cv_data$`# non-zero Loadings`[which.min(cv_data$`CV Errors`)]
    best_ridge <- cv_data$`Ridge Penalty`[which.min(cv_data$`CV Errors`)]
    best_cv_error <- min(cv_data$`CV Errors`)
    y_min_limit <- best_cv_error
    y_max_limit <- max(cv_data$`CV Errors`)
    out_list$`CV Results` <- ggplot(cv_data, aes(x = `Ridge Penalty`, y = `CV Errors`, colour = `Sparsity Ratio`, 
                                                 shape = `Sparsity Ratio`)) +
      geom_point(size = 2) +
      geom_hline(yintercept = cv_data$`CV Errors`[1], colour = "black") +
      scale_colour_manual(values =  c("#88CCEE", "#44799E", "#000000", "#41784A", "#117733"),
                          name = "Sparsity Ratio") +
      scale_shape_discrete(name = "Sparsity Ratio") +
      geom_point(data = subset(cv_data, `CV Errors` == min(`CV Errors`)), aes(x = `Ridge Penalty`, y = `CV Errors`), 
                 colour = "black", fill = "#882255", size = 3, shape = 22) +
      scale_y_continuous(trans = "log10", limits = c(y_min_limit, y_max_limit)) + 
      scale_x_continuous(trans = "log10") +
      annotate("text",  x = best_ridge, y = best_cv_error,
               label = best_combo, angle = 0, vjust = 1.6, hjust = 1, size = 4, color = "darkred") +
      labs(x = "log Ridge Penalty",
           y = "log CV Error") +
      theme_minimal()
    
    # BIC results
    bic_data <- data.frame(x$BIC$`BIC Results`, check.names = FALSE)
    bic_data$`# non-zero Loadings` <- c(NaN)
    for(i in 1:dim(bic_data)[1]){
      bic_data$`# non-zero Loadings`[i] <- paste0("(", 
                                                  paste0(sprintf("%.0f", bic_data[i, 3:(dim(bic_data)[2] - 1)]), collapse = ";"), 
                                                  ")")
    }
    bic_data$`Sparsity Ratio` <- as.factor(binned_data_equal_width)
    best_bic_combo <- bic_data$`# non-zero Loadings`[which.min(bic_data$`BIC`)]
    best_bic_ridge <- bic_data$`Ridge Penalty`[which.min(bic_data$`BIC`)]
    best_bic <- min(bic_data$`BIC`)
    y_bic_min_limit <- best_cv_error
    y_bic_max_limit <- max(bic_data$`BIC`)
    out_list$`BIC Results` <- ggplot(bic_data, aes(x = `Ridge Penalty`, y = `BIC`, colour = `Sparsity Ratio`, 
                                                   shape = `Sparsity Ratio`)) +
      geom_point(size = 2) +
      geom_hline(yintercept = bic_data$BIC[1], colour = "black") +
      scale_colour_manual(values =  c("#88CCEE", "#44799E", "#000000", "#41784A", "#117733"),
                          name = "Sparsity Ratio") +
      scale_shape_discrete(name = "Sparsity Ratio") +
      geom_point(data = subset(bic_data, `BIC` == min(`BIC`)), aes(x = `Ridge Penalty`, y = `BIC`), 
                 colour = "black", fill = "#882255", size = 3, shape = 22) +
      scale_y_continuous(limits = c(y_bic_min_limit, y_bic_max_limit)) + 
      scale_x_continuous(trans = "log10") +
      annotate("text",  x = best_bic_ridge, y = best_bic,
               label = best_bic_combo, angle = 0, vjust = 1.6, hjust = 1, size = 4, color = "darkred") +
      labs(x = "log Ridge Penalty",
           y = "BIC") +
      theme_minimal()
    
  }else if(any(grepl("Maximum No. of LARS Steps", colnames(x$CV$`CV Results`), ignore.case = FALSE))){
    cv_data <- data.frame(x$CV$`CV Results`, check.names = FALSE)
    breaks <- floor(seq(from = min(cv_data$`Maximum No. of LARS Steps`, na.rm = TRUE), 
                        to = max(cv_data$`Maximum No. of LARS Steps`, na.rm = TRUE), 
                        length.out = 6))
    breaks[length(breaks)] <- breaks[length(breaks)] + 0.000006
    labels <- paste0("[", sprintf("%.0f", floor(breaks[1:5] / 0.01) * 0.01),
                     ", ",
                     sprintf("%.0f", floor(breaks[2:6] / 0.01) * 0.01),
                     ")")
    labels[5] <- paste0("[", sprintf("%.0f", floor(breaks[5] / 0.01) * 0.01),
                        ", ",
                        sprintf("%.0f", floor(breaks[6] / 0.01) * 0.01),
                        "]")
    binned_data_equal_width <- cut(cv_data$`Maximum No. of LARS Steps`,
                                   breaks = breaks,
                                   right = FALSE,
                                   include.lowest = TRUE,
                                   labels = labels)
    cv_data$`# of LARS Steps` <- as.factor(binned_data_equal_width)
    best_combo <- cv_data$`Maximum No. of LARS Steps`[which.min(cv_data$`CV Errors`)]
    best_ridge <- cv_data$`Ridge Penalty`[which.min(cv_data$`CV Errors`)]
    best_cv_error <- min(cv_data$`CV Errors`)
    y_min_limit <- best_cv_error
    y_max_limit <- max(cv_data$`CV Errors`)
    out_list$`CV Results` <- ggplot(cv_data, aes(x = `Ridge Penalty`, y = `CV Errors`, colour = `# of LARS Steps`, 
                                                 shape = `# of LARS Steps`)) +
      geom_point(size = 2) +
      geom_hline(yintercept = cv_data$`CV Errors`[1], colour = "black") +
      scale_colour_manual(values =  c("#88CCEE", "#44799E", "#000000", "#41784A", "#117733"),
                          name = "# of LARS Steps") +
      scale_shape_discrete(name = "# of LARS Steps") +
      geom_point(data = subset(cv_data, `CV Errors` == min(`CV Errors`)), aes(x = `Ridge Penalty`, y = `CV Errors`), 
                 colour = "black", fill = "#882255", size = 3, shape = 22) +
      scale_y_continuous(trans = "log10", limits = c(y_min_limit, y_max_limit)) + 
      scale_x_continuous(trans = "log10") +
      annotate("text",  x = best_ridge, y = best_cv_error,
               label = best_combo, angle = 0, vjust = 1.6, hjust = 1, size = 4, color = "darkred") +
      labs(x = "log Ridge Penalty",
           y = "log CV Error") +
      theme_minimal()
    
    # BIC results
    bic_data <- data.frame(x$BIC$`BIC Results`, check.names = FALSE)
    bic_data$`# of LARS Steps` <- as.factor(binned_data_equal_width)
    best_bic_combo <- bic_data$`Maximum No. of LARS Steps`[which.min(bic_data$`BIC`)]
    best_bic_ridge <- bic_data$`Ridge Penalty`[which.min(bic_data$`BIC`)]
    best_bic <- min(bic_data$`BIC`)
    y_bic_min_limit <- best_cv_error
    y_bic_max_limit <- max(bic_data$`BIC`)
    out_list$`BIC Results` <- ggplot(bic_data, aes(x = `Ridge Penalty`, y = `BIC`, colour = `# of LARS Steps`, 
                                                   shape = `# of LARS Steps`)) +
      geom_point(size = 2) +
      geom_hline(yintercept = bic_data$BIC[1], colour = "black") +
      scale_colour_manual(values =  c("#88CCEE", "#44799E", "#000000", "#41784A", "#117733"),
                          name = "# of LARS Steps") +
      scale_shape_discrete(name = "# of LARS Steps") +
      geom_point(data = subset(bic_data, `BIC` == min(`BIC`)), aes(x = `Ridge Penalty`, y = `BIC`), 
                 colour = "black", fill = "#882255", size = 3, shape = 22) +
      scale_y_continuous(limits = c(y_bic_min_limit, y_bic_max_limit)) + 
      scale_x_continuous(trans = "log10") +
      annotate("text",  x = best_bic_ridge, y = best_bic,
               label = best_bic_combo, angle = 0, vjust = 1.6, hjust = 1, size = 4, color = "darkred") +
      labs(x = "log Ridge Penalty",
           y = "BIC") +
      theme_minimal()
  }
  
  return(out_list)
}
