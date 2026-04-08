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

#' @name twoStepSDFM
#' @title Two Step Sparse Dynamic Factor Model Estimator.
#' @description
#' Estimate a sparse dynamic factor model with measurement equation
#' \deqn{
#'   \bm{x}_t = \bm{\Lambda} \bm{f}_{t} + \bm{\xi}_t,\quad \bm{\xi}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{\xi}),
#' }
#' and transition equation
#' \deqn{
#'   \bm{f}_t = \sum_{p=0}^P\bm{\Phi}_p \bm{f}_{t-p} + \bm{\epsilon}_t,\quad \bm{\epsilon}_t \sim \mathcal{N}(\bm{0}, \bm{\Sigma}_{f}).
#' }
#' using sparse principal components  analysis and the Kalman Filter and 
#' Smoother according to \insertRef{franjic2024nowcasting}{TwoStepSDFM}.
#' 
#' @param data Numeric (no_of_vars \eqn{\times}{x} no_of_obs) matrix of data or 
#' zoo/xts object sampled at the same frequency.
#' @param delay Integer vector of variable delays, measured as the number of 
#' months since the latest available observation.
#' @param selected Integer vector of the number of selected variables for each 
#' factor.
#' @param no_of_factors Integer number of factors.
#' @param max_factor_lag_order Integer maximum order of the VAR process in the 
#' transition equation.
#' @param lag_estim_criterion Information criterion used for the estimation of 
#' the factor VAR order (`"BIC"` (default), `"AIC"`, `"HIC"`).
#' @param decorr_errors Logical, whether or not the errors should be 
#' decorrelated.
#' @param ridge_penalty Ridge penalty.
#' @param lasso_penalty Numeric vector, lasso penalties for each factor (set to 
#' NULL to disable as stopping criterion).
#' @param max_iterations Integer maximum number of iterations.
#' @param max_no_steps Integer number of LARS steps (set to NULL to disable as 
#' stopping criterion).
#' @param weights Numeric vector, weights for each variable weighing the 
#' \eqn{\ell_1}{`l_1`} size constraint.
#' @param comp_null Numeric computational zero.
#' @param spca_conv_crit Conversion threshold for the SPCA algorithm.
#' @param parallel Logical, whether or not to use Eigen's internal parallel 
#' matrix operations.
#' @param fcast_horizon Integer number of additional Filter predictions into the 
#' future.
#' @param jitter Numerical jitter for stability of internal solver algorithms. 
#' The jitter is added to the diagonal entries of the variance covariance matrix 
#' of the measurement errors.
#' @param svd_method Either "fast" or "precise". Option "fast" uses Eigen's 
#' BDCSVD divide and conquer method for the computation of the singular values. 
#' Option "precise" (default) implements the slower, but numerically more stable 
#' JacobiSVD method \insertCite{eigenweb}{TwoStepSDFM}.
#' 
#' @details
#' The function performs a two-step estimation procedure for sparse dynamic 
#' factor models as described in \insertRef{franjic2024nowcasting}{TwoStepSDFM}. 
#' In the first step, the factor loading matrix is estimated using SPCA 
#' \insertCite{zou2006sparse}{TwoStepSDFM}. This will shrink some of the 
#' loadings towards or exactly to zero. In the second step the latent factors
#' are estimated using the univariate representation of the Kalman Filter and
#' Smoother \insertCite{koopman2000fast}{TwoStepSDFM}.
#' 
#' The function takes three stopping criteria for the SPCA algorithm:
#' `selected`, `lasso_penalty`, and `max_no_steps`. The argument `weights`
#' allows specifying weights for the \eqn{\ell_1}{l1} constraint.  `svd_method` 
#' controls the decomposition method for internal SVDs. For a detailed 
#' description of these arguments and the SPCA step, see 
#' \code{\link{sparsePCA}}.
#' 
#' With respect to the univariate representation of the Kalman filter and 
#' smoother, `decorr_errors` indicates whether the data should be decorrelated 
#' internally prior to filtering and smoothing. `jitter` is added to the 
#' diagonal elements of the measurement variance–covariance matrix. For more 
#' details, see \code{\link{kalmanFilterSmoother}}.
#' 
#' For more information on the two-step estimation procedure see 
#' \insertRef{franjic2024nowcasting}{TwoStepSDFM}.
#' 
#' @return 
#' An object of class SDFMFit with main components:
#' #' \describe{
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
#' \insertRef{zou2006sparse}{TwoStepSDFM}
#' 
#' \insertRef{eigenweb}{TwoStepSDFM}
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
#' \code{\link{twoStepDenseDFM}}: Two-step estimation routine for a dense 
#' dynamic factor model.
#' 
#' @examples
#' data(factor_model)
#' no_of_vars <- dim(factor_model$data)[2]
#' no_of_factors <- dim(factor_model$factors)[2]
#' sdfm_fit <- twoStepSDFM(data = factor_model$data, delay = factor_model$delay,
#'                         selected = rep(floor(0.5 * no_of_vars), no_of_factors),
#'                         no_of_factors = no_of_factors)
#' print(sdfm_fit)
#' sdfm_plots <- plot(sdfm_fit)
#' sdfm_plots$`Factor Time Series Plots`
#' sdfm_plots$`Loading Matrix Heatmap`
#' sdfm_plots$`Meas. Error Var.-Cov. Matrix Heatmap`
#' sdfm_plots$`Meas. Error Var.-Cov. Eigenvalue Plot`
#' 
#' @export
twoStepSDFM <- function (data, 
                         delay, 
                         selected, 
                         no_of_factors, 
                         max_factor_lag_order = 10, 
                         lag_estim_criterion = "BIC", 
                         decorr_errors = TRUE, 
                         ridge_penalty = 1e-06, 
                         lasso_penalty = NULL, 
                         max_iterations = 1000, 
                         max_no_steps = NULL, 
                         weights = NULL, 
                         comp_null = 1e-15, 
                         spca_conv_crit = 1e-04, 
                         parallel = FALSE, 
                         fcast_horizon = 0, 
                         jitter = 1e-08, 
                         svd_method = "precise") {
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
    stop(paste0("Too few observations as no_of-variables >= no_of_observations."))
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
    stop(paste0("max_factor_lag_order must be strictly positve."))
  }
  decorr_errors <- checkBoolean(decorr_errors, "decorr_errors")
  if (is.null(lag_estim_criterion)) {
    stop(paste0("lag_estim_criterion must be either \"BIC\", \"AIC\", or \"HIC\"."))
  }
  if (!(lag_estim_criterion %in% c("AIC", "BIC", "HIC"))) {
    stop(paste0("lag_estim_criterion must be either \"BIC\", \"AIC\", or \"HIC\"."))
  }
  spca_conv_crit <- checkPositiveDouble(spca_conv_crit, "spca_conv_crit")
  if (spca_conv_crit == 0) {
    warning("spca_conv_crit should not be exactly 0. It will be jittered before further use.")
    spca_conv_crit <- 1e-15
  }
  parallel <- checkBoolean(parallel, "parallel")
  fcast_horizon <- checkPositiveSignedInteger(fcast_horizon, "fcast_horizon")
  jitter <- checkPositiveDouble(jitter, "jitter")
  if (!(svd_method %in% c("fast", "precise"))) {
    stop(paste0("svd_method must be \"fast\" for usage of Eigen's BDCSVD or \"precise\" for usage of Eigen's JacobiSVD"))
  }
  comp_null <- checkPositiveDouble(comp_null, "comp_null")
  if (comp_null == 0) {
    warning("comp_null should not be exactly 0. It will be jittered before further use.")
    comp_null <- 1e-15
  }
  
  # Mishandling of selected
  if (is.null(selected)) {
    selected <- matrix(rep(no_of_variables, no_of_factors), ncol = 1)
  }
  else {
    selected <- checkPositiveSignedParameterVector(selected, "selected", no_of_factors)
  }
  if (any(selected > no_of_variables)) {
    warning(paste0("The elements in selected should not exceed the number of variables ", 
                   no_of_variables, ". The corresponding variables are set to ", no_of_variables, "."))
    selected[which(selected > no_of_variables), ] <- no_of_variables
  }
  
  # Mishandling of ridge_penalty and lasso_penalty
  ridge_penalty <- checkPositiveDouble(ridge_penalty, "ridge_penalty")
  if (!is.null(lasso_penalty)) {
    if (!is.numeric(lasso_penalty) || any(is.na(lasso_penalty))) {
      stop(paste0("lasso_penalty must be a vector of non-NA numeric values."))
    }
    if (length(lasso_penalty) != no_of_factors) {
      stop(paste0("lasso_penalty must be of length no_of_factors = ", no_of_factors))
    }
    if (any(lasso_penalty < 0)) {
      stop(paste0("All elements of lasso_penalty non-negative."))
    }
  }
  else {
    lasso_penalty <- rep(-2147483647L, no_of_factors)
  }
  
  # Mishandling of max_iterations and max_no_steps
  max_iterations <- checkPositiveSignedInteger(max_iterations, "max_iterations")
  if (!is.null(max_no_steps)) {
    max_no_steps <- checkPositiveSignedInteger(max_no_steps, "max_no_steps")
    if (max_no_steps == 0) {
      stop(paste0("max_no_steps must be strictly positve."))
    }
  }
  else {
    max_no_steps <- -2147483647L
  }
  
  # Mishandling of weights
  if (!is.null(weights)) {
    if (!is.numeric(weights) || any(is.na(weights)) || any(weights <= 
                                                           0)) {
      stop(paste0("weights must be a vector of non-NA numeric values strictly greater 0."))
    }
    if (length(weights) != no_of_variables) {
      stop(paste0("weights must be of length no_of_variables = ", no_of_variables))
    }
    if (sum(weights) != 1) {
      message("weights are standardised to sum to 1")
      weights <- weights / sum(weights)
    }
  }
  else {
    weights <- rep(1, no_of_variables)
  }
  
  result <- runSDFMKFS(X_in = data_r, delay = delay, selected = selected, 
                       R = as.integer(no_of_factors), order = as.integer(max_factor_lag_order), 
                       decorr_errors = decorr_errors, crit = lag_estim_criterion, 
                       l2 = ridge_penalty, l1 = lasso_penalty, max_iterations = as.integer(max_iterations), 
                       steps = max_no_steps, weights = weights, comp_null = comp_null, 
                       spca_conv_crit = spca_conv_crit, parallel = parallel, 
                       fcast_horizon = fcast_horizon, jitter = jitter, svd_method = svd_method)
  
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

#' @name print.SDFMFit
#' @title Generic printing function for SDFMFit S3 objects
#' @description
#' Print a compact summary of an `SDFMFit` object.
#'
#' @param x `SDFMFit` object.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' No return value; Prints a summary to the console.
#'
#' @author
#' Domenic Franjic
#' 
#' @export
print.SDFMFit <- function (x, ...) 
{
  simulated_time_series <- is.zoo(x$smoothed_factors)
  no_of_factors <- ifelse(simulated_time_series, dim(x$smoothed_factors)[2], dim(x$smoothed_factors)[1])
  no_of_obs <- ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2])
  cat("Simulated Dynamic Factor Model\n")
  cat("=========================================================================\n")
  cat("No. of Observations                        :", ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2]), "\n")
  cat("No. of Variables                           :", ifelse(simulated_time_series, dim(x$data)[2], dim(x$data)[1]), "\n")
  cat("No. of Factors                             :", no_of_factors, "\n")
  cat("Factor Lag Order                           :", x$factor_var_lag_order,  "\n")
  cat("No. of zero elements in the loading matrix :", sum(x$loading_matrix_estimate == 0), "\n")
  if (x$llt_success_code == -1) {
    cat("Info: LLT failed. Used robust LDLT instead.\n")
  }
  else if (x$llt_success_code == -1) {
    cat("Warning: LLT and LDLT failed. Used uncorrelated errors.\n")
  }
  cat("=========================================================================\n")
  cat("Head of the factors :\n")
  max_print <- min(5, no_of_obs)
  if (simulated_time_series) {
    print(head(x$smoothed_factors, max_print))
  }
  else {
    print(x$smoothed_factors[, 1:max_print])
  }
  cat("Tail of the factors :\n")
  if (simulated_time_series) {
    print(tail(x$smoothed_factors, max_print))
  }
  else {
    print(x$smoothed_factors[, (dim(x$smoothed_factors)[2] - (max_print - 1)):(dim(x$smoothed_factors)[2])])
  }
  max_print_loadings <- min(5, ifelse(simulated_time_series, dim(x$smoothed_factors)[1], dim(x$smoothed_factors)[2]))
  cat("Head of the loading matrix :\n")
  print(head(x$loading_matrix_estimate, max_print_loadings))
  cat("Tail of the loading matrix :\n")
  print(tail(x$loading_matrix_estimate, max_print_loadings))
  cat("=========================================================================\n")
  invisible(x)
}

#' @name plot.SDFMFit
#' @title Generic plotting function for SDFMFit S3 objects
#' @param x `SDFMFit` object.
#' @param axis_text_size Numeric size of x- and y-axis labels. Passed to ggplot2 
#' `theme(..., text = element_text(size = axis_text_size))`.
#' @param legend_title_text_size Numeric size of x- and y-axis labels. Passed to
#' ggplot2 
#' `theme(..., legend.title = element_text(size = legend_title_text_size))`.
#' @param ... Additional parameters for the plotting functions.
#' 
#' @return
#' A named list of plot objects:
#' \describe{
#'   \item{`Factor Time Series Plots`}{`patchwork`/`ggplot` object graphing the 
#'   estimated factors over time with 95% confidence bands based on the smoother 
#'   uncertainty of the Kalman Filter and Smoother.}
#'   \item{`Loading Matrix Heatmap`}{`ggplot` object showing a heatmap of the 
#'   estimated factor loadings. Zeros are highlighted in black.}
#'   \item{`Meas. Error Var.-Cov. Matrix Heatmap`}{`ggplot` object showing a 
#'   heatmap of the measurement error variance-covariance matrix.}
#'   \item{`Eigenvalue Plot`}{`ggplot` object showing a bar plot of the 
#'   eigenvalues of the measurement error variance–covariance matrix.}
#'  }
#' 
#' @author
#' Domenic Franjic
#' 
#' @export
plot.SDFMFit <- function (x, 
                          axis_text_size = 20, 
                          legend_title_text_size = 20, 
                          ...) 
{
  out_list <- list()
  if (is.zoo(x$data)) {
    series_names <- colnames(x$data)
    no_of_factors <- dim(x$smoothed_factors)[2]
    no_of_obs <- dim(x$factor_estimate)[2]
    time_vector <- as.Date(time(x$smoothed_factors))
    factors <- x$smoothed_factors
  }
  else {
    series_names <- rownames(x$data)
    no_of_factors <- dim(x$smoothed_factors)[1]
    no_of_obs <- dim(x$factor_estimate)[2]
    time_vector <- 1:dim(x$smoothed_factors)[2]
    factors <- t(x$smoothed_factors)
    factors <- as.zoo(ts(factors, start = c(1, 1), frequency = 12))
  }
  
  out_list$`Factor Time Series Plots` <- plotFactorEstimates(factors, x$smoothed_state_variance, no_of_factors, axis_text_size)
  
  out_list$`Loading Matrix Heatmap` <- plotLoadingHeatMap(x$loading_matrix_estim, 
                                                          series_names, no_of_factors,
                                                          axis_text_size, legend_title_text_size)
  
  if (is.character(x$error_var_cov_cholesky_factor)) {
    if (is.zoo(x$data)) {
      residuals <- coredata(na.omit(x$data)) - coredata(x$factor_estimate) %*% t(x$loading_matrix_estim)
    }
    else {
      residuals <- na.omit(t(x$data)) - t(x$factor_estimate) %*% t(x$loading_matrix_estim)
    }
    measurement_error_var_cov_df <- as.data.frame(t(residuals) %*% residuals * 1/(dim(residuals)[1] - 1))
    out_list$`Meas. Error Var.-Cov. Matrix Heatmap` <- plotMeasVarCovHeatmap(measurement_error_var_cov_df, 
                                                                             series_names, 
                                                                             axis_text_size, 
                                                                             legend_title_text_size)
  }
  else {
    measurement_error_var_cov_df <- as.data.frame(x$error_var_cov_cholesky_factor %*% t(x$error_var_cov_cholesky_factor))
    out_list$`Meas. Error Var.-Cov. Matrix Heatmap` <- plotMeasVarCovHeatmap(measurement_error_var_cov_df, 
                                                                             series_names, 
                                                                             axis_text_size, 
                                                                             legend_title_text_size)
  }
  
  out_list$`Meas. Error Var.-Cov. Eigenvalue Plot` <- plotMeasVarCovEigenvalues(eigen(measurement_error_var_cov_df)$values, 
                                                                                no_of_factors, axis_text_size, legend_title_text_size)
  
  return(out_list)
}

#' @name predict.SDFMFit
#' @title Generic plotting function for SDFMFit S3 objects
#' @description
#' Predict all missing observations due to ragged edges in the data set plus 
#' horizon steps ahead.
#' 
#' @param object `SDFMFit` object.
#' @param horizon Number of forecasting steps into the future. Must be smaller 
#' than or equal to `x$factor_fcast_horizon`.
#' @param ... Additional parameters for the prediction function.
#' 
#' @return
#' A named list of plot objects:
#' \describe{
#'   \item{data}{Object containing the original data. The object inherits its 
#'   class from `object$data`: If data is provided as `zoo`, `data` will be a 
#'   `zoo` object. If `data` is provided as matrix, `data` will be a 
#'   (`no_of_factors`\eqn{\times}{x}`no_of_obs`) matrix.}
#'   \item{data_missing_pred}{Object containing only the predictions of all 
#'   missing observations plus the forecasts. Inherits its class from 
#'   `object$data` as above.}
#'   \item{data_imputed}{Object containing the observed data, predictions of 
#'   all missing observations plus the forecasts. Inherits its class from 
#'   `object$data` as above.}
#'  }
#' 
#' @author
#' Domenic Franjic
#' 
#' @examples
#' data(factor_model)
#' no_of_vars <- dim(factor_model$data)[2]
#' no_of_factors <- dim(factor_model$factors)[2]
#' sdfm_fit <- twoStepSDFM(data = factor_model$data, delay = factor_model$delay,
#'                         selected = rep(floor(0.5 * no_of_vars), no_of_factors),
#'                         no_of_factors = no_of_factors, fcast_horizon = 5)
#' dfm_fit <- twoStepDenseDFM(data = factor_model$data, delay = factor_model$delay, 
#'                            no_of_factors = no_of_factors, fcast_horizon = 5)
#' predict(sdfm_fit, horizon = 5)
#' predict(dfm_fit, horizon = 5)
#' 
#' @export
predict.SDFMFit <- function (object, 
                             horizon = 0,
                             ...) {
  horizon <- checkPositiveSignedInteger(horizon, "horizon")
  if (horizon > object$factor_fcast_horizon) {
    stop("There are not enough forecasts of the estimated factors. Re-run twoStepSDFM/twoStepDenseDFM setting fcast_horizon >= horizon")
  }
  if (is.zoo(object$data)) {
    no_of_obs <- dim(object$data)[1]
    no_of_vars <- dim(object$data)[2]
    time_vector <- as.Date(time(object$smoothed_factors))
    factors <- object$smoothed_factors
    data <- object$data
  }
  else {
    no_of_obs <- dim(object$data)[2]
    no_of_vars <- dim(object$data)[1]
    time_vector <- 1:dim(object$smoothed_factors)[2]
    factors <- t(object$smoothed_factors)
    factors <- as.zoo(ts(factors, start = c(1, 1), frequency = 12))
    data <- as.zoo(ts(data, start = c(1, 1), frequency = 12))
  }
  
  data_pred <- (object$loading_matrix_estimate %*% t(coredata(factors)))[, 1:(no_of_obs + horizon)]
  object$data
  result <- list()
  result$data <- object$data
  result$data_missing_pred <- matrix(NaN, no_of_vars, no_of_obs + horizon)
  result$data_imputed <- matrix(NaN, no_of_vars, no_of_obs + horizon)
  result$data_imputed[, 1:no_of_obs] <- t(coredata(data))
  for (var in 1:no_of_vars) {
    if (object$data_delay[var] > 0) {
      curr_predictions <- data_pred[var, (no_of_obs - object$data_delay[var] + 1):(no_of_obs + horizon)]
      result$data_missing_pred[var, (no_of_obs - object$data_delay[var] + 1):(no_of_obs + horizon)] <- curr_predictions
      result$data_imputed[var, (no_of_obs - object$data_delay[var] + 1):(no_of_obs + horizon)] <- curr_predictions
    }
    else if (horizon > 0) {
      curr_predictions <- data_pred[var, (no_of_obs + 1):(no_of_obs + horizon)]
      result$data_missing_pred[var, (no_of_obs + 1):(no_of_obs + horizon)] <- curr_predictions
      result$data_imputed[var, (no_of_obs + 1):(no_of_obs + horizon)] <- curr_predictions
    }
  }
  if (is.zoo(object$data)) {
    result$data_missing_pred <- as.zoo(ts(t(result$data_missing_pred), start = c(year(time_vector[1]), month(time_vector[1])), 
                                          frequency = 12))
    result$data_imputed <- as.zoo(ts(t(result$data_imputed), start = c(year(time_vector[1]), month(time_vector[1])), 
                                     frequency = 12))
  }
  
  return(result)
}