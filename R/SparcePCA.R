#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @import zoo
#' @import xts
#' @import lubridate
#' @import ggplot2
#' @import stats
#' @import utils
#' @import grDevices
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

#' @name sparsePCA
#' @title Sparse Principal Components Analysis
#' @description
#' Estimate sparse sparse principal components via SPCA according to 
#' \insertRef{zou2006sparse}{TwoStepSDFM}.
#' 
#' @param data Numeric (no_of_vars \eqn{\times}{x} no_of_obs) matrix of data or 
#' zoo/xts object sampled at the same frequency.
#' @param delay Integer vector of variable delays, measured as the number of 
#' months since the latest available observation.
#' @param selected Integer vector of the number of selected variables for each 
#' factor.
#' @param no_of_factors Integer number of factors.
#' @param ridge_penalty Numeric ridge penalty.
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
#' @param svd_method Either "fast" or "precise". Option "fast" uses Eigen's 
#' BDCSVD divide and conquer method for the computation of the singular values. 
#' Option "precise" (default) implements the slower, but numerically more stable 
#' JacobiSVD method.
#' @param normalise Logical, whether to normalise the loading matrix as in 
#' \insertRef{zou2020elnet}{TwoStepSDFM}. Default is `TRUE`.
#' @param comp_var_expl Logical, whether to compute the relative variance 
#' explained by each factor. Default is `TRUE``.`
#'
#' @details
#' The function takes three stopping criteria: `selected`, `lasso_penalty`, and
#' `max_no_steps`. With `selected` the SPCA algorithm stops if each column of 
#' the estimated loading matrix has the corresponding number of non-zero 
#' loadings. This allows the user to directly control the degree of sparsity of 
#' each factor loading. With `lasso_penalty`, the SPCA algorithm stops as soon 
#' as the side-constraints of the inherent elastic-net problem are no longer 
#' satisfied. With `max_no_steps`, the SPCA algorithm only takes that many LARS 
#' steps for each factor loading's individual elastic-net problem before 
#' stopping. If all criteria are provided, the first one satisfied will stop the
#' algorithm. For details see also \insertCite{zou2006sparse}{TwoStepSDFM} and 
#' \insertCite{zou2020elnet}{TwoStepSDFM}.
#' 
#' Loosely, each SPCA algorithm iteration solves an elastic-net type problem for
#' each column of the loading matrix. One can extend this problem to the 
#' adaptive elastic-net \insertCite{zou2009adaptive}{TwoStepSDFM}. The variable 
#' `weights` lets the user provide weights for each observation. These weights 
#' must be strictly greater than zero and are normalised internally to represent 
#' relative weights. For more information on the computational implementation of
#' the weight extension in the context of SPCA see 
#' \insertRef{zou2024general}{TwoStepSDFM}.
#'  
#' In each SPCA algorithm iteration, the function executes an SVD. To this end, 
#' Eigen provides two alternatives \insertCite{eigenweb}{TwoStepSDFM}: Option 
#' `precise` makes use of JacobiSVD. This method is numerically more stable, but 
#' computationally costly, especially for medium to large matrices. Option 
#' `fast` makes use of BDCSVD. This divide-and-conquer approach can lead to 
#' significant performance gains with respect to large matrices. BDCSVD, 
#' however, can be numerically unstable when Eigen is compiled with aggressive 
#' speed optimisations. In the context of the `R`, this should be of no concern.
#' By default, `R` and most packages are compiled with "mild" `-O2` optimisation
#' and without any additional aggressive optimisation flags. Nonetheless, one 
#' should checker whether both variants provide reasonably close results before 
#' switching to `fast`. For more information see 
#' \insertRef{eigenweb}{TwoStepSDFM}.
#' 
#' @return
#' An object of class `SPCAFit` with components:
#' \describe{
#'   \item{data}{Original data matrix.}
#'   \item{loading_matrix_estim}{Numeric matrix of estimated factor loadings.}
#'   \item{factor_estim}{Object containing the SPCA factor estimates. The 
#'   object inherits its class from `data`: If `data` is provided as `zoo`, 
#'   `factor_estim` will be a `zoo` object. If `data` is provided as  `matrix`, 
#'   `factor_estim` will be a (`no_of_factors` \eqn{\times}{x} `no_of_obs`) 
#'   matrix.}
#'   \item{total_var_expl}{Numeric total variance explained.}
#'   \item{pct_var_expl}{Numeric vector relative variance explained by each 
#'   factor.}
#' }
#'
#' @author
#' Domenic Franjic
#' 
#' @references
#' \insertRef{zou2006sparse}{TwoStepSDFM}
#' 
#' \insertRef{zou2009adaptive}{TwoStepSDFM}
#' 
#' \insertRef{eigenweb}{TwoStepSDFM}
#' 
#' \insertRef{zou2020elnet}{TwoStepSDFM}
#' 
#' \insertRef{zou2024general}{TwoStepSDFM}
#'
#' @examples
#' data(factor_model)
#' set.seed(17032026)
#' no_of_factors <- 3
#' no_of_vars <- dim(factor_model$data)[2]
#' selected <- rep(floor(0.5 * no_of_vars), no_of_factors)
#' lasso_penalty <- exp(runif(no_of_factors, -10, 1))
#' max_no_steps <- 1000
#' spca_fit <- sparsePCA(data = factor_model$data, delay = factor_model$delay, selected = selected,
#'                       no_of_factors = no_of_factors, ridge_penalty = 1e-2, lasso_penalty = lasso_penalty,
#'                       max_iterations = 1000, weights = NULL, max_no_steps = max_no_steps, comp_null = 1e-15,
#'                       spca_conv_crit = 1e-04, parallel = FALSE, svd_method = "precise", normalise = FALSE,
#'                       comp_var_expl = TRUE
#' )
#' print(spca_fit)
#' spca_plots <- plot(spca_fit)
#' spca_plots$`Factor Time Series Plots`
#' spca_plots$`Loading Matrix Heatmap`
#' spca_plots$`Meas. Error Var.-Cov. Matrix Heatmap`
#' spca_plots$`Eigenvalue Plot`
#' spca_plots$`Variance Explained Chart`
#'
#' @export
sparsePCA <- function(data,
                      delay,
                      selected,
                      no_of_factors,
                      ridge_penalty = 1e-6,
                      lasso_penalty = NULL,
                      max_iterations = 1000,
                      weights = NULL,
                      max_no_steps = NULL,
                      comp_null = 1e-15,
                      spca_conv_crit = 1e-4,
                      parallel = FALSE,
                      svd_method = "precise",
                      normalise = TRUE,
                      comp_var_expl = TRUE) {
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
  no_of_variables <- dim(data_r)[2]
  no_of_observations <- dim(data_r)[1]
  if(is.null(delay)){
    delay <- matrix(rep(0, no_of_variables), ncol = 1)
  }else{
    delay <- checkPositiveSignedParameterVector(delay, "delay", no_of_variables)
  }
  
  # Check for NAs in the dataset outside the ragged edges
  na_ind <- FALSE
  for(col in 1:dim(data_r)[2]){
    na_ind <- any(is.na(data_r[1:(no_of_observations - delay[col]), col]))
  }
  if(na_ind){
    warning(paste0("data has NA values. AccordingThe corresponding time points will not be considered in th estimation of the loading matrix.")) 
  }
  
  # # Misshandling of dimensions
  # if(no_of_variables >= no_of_observations){
  #   warning(paste0("Too few observations as no_of-variables >= no_of_observations."))
  # }
  
  # Mishandling of selected
  if(is.null(selected)){
    selected <- matrix(rep(no_of_variables, no_of_factors), ncol = 1)
  }else{
    selected <- checkPositiveSignedParameterVector(selected, "selected", no_of_factors)
  }
  if(any(selected > no_of_variables)){
    warning(paste0("The elements in selected should not exceed the number of variables ", no_of_variables, ". The corresponding variables are set to ", no_of_variables, "."))
    selected[which(selected > no_of_variables),] <- no_of_variables
  }
  
  # Mishandling of number of factors
  no_of_factors <- checkPositiveSignedInteger(no_of_factors, "no_of_factors")
  if(no_of_factors == 0){
    stop("no_of_factors must be strictly positive.")
  }
  if(no_of_factors > no_of_variables){
    stop(paste0("no_of_factors must be smaller than no_of_variables."))
  }
  
  # Mishandling of ridge penalty
  ridge_penalty <- checkPositiveDouble(ridge_penalty, "ridge_penalty")
  
  # Mishandling of lasso_penalty penalty
  if (!is.null(lasso_penalty)){
    if(!is.numeric(lasso_penalty) || any(is.na(lasso_penalty))){
      stop(paste0("lasso_penalty must be a vector of non-NA numeric values."))
    }
    if(length(lasso_penalty) != no_of_factors){
      stop(paste0("lasso_penalty must be of length no_of_factors = ", no_of_factors))
    }
    if(any(lasso_penalty < 0)){
      stop(paste0("All elements of lasso_penalty non-negative."))
    }
  }else{
    lasso_penalty <- rep(-2147483647L, no_of_factors)
  }
  
  # Mishandlilng of max_iterations
  max_iterations <- checkPositiveSignedInteger(max_iterations, "max_iterations")
  
  # Mishandling of max_no_steps
  if(!is.null(max_no_steps)){
    max_no_steps <- checkPositiveSignedInteger(max_no_steps, "max_no_steps")
    if(max_no_steps == 0){
      stop(paste0("max_no_steps must be strictly positve."))
    }
  }else{
    max_no_steps <- -2147483647L # C++ INT_MIN
  }
  
  if(!is.null(weights)){
    if(!is.numeric(weights) || any(is.na(weights)) || any(weights <= 0)){
      stop(paste0("weights must be a vector of non-NA numeric values strictly greater 0."))
    }
    if(length(weights) != no_of_variables){
      stop(paste0("weights must be of length no_of_variables = ", no_of_variables))
    }
    if(sum(weights) != 1){
      message("weights are standardised to sum to 1")
      weights <- weights / sum(weights)
    }
  }else{
    weights <- rep(1, no_of_variables)
  }
  
  # Mishandling of comp_null
  comp_null <- checkPositiveDouble(comp_null, "comp_null")
  if(comp_null == 0){
    warning("comp_null should not be exactly 0. It will be jittered before further use.")
    comp_null <- 1e-15
  }
  
  # Mishandling of spca_conv_crit
  spca_conv_crit <- checkPositiveDouble(spca_conv_crit, "spca_conv_crit")
  if(spca_conv_crit == 0){
    warning("spca_conv_crit should not be exactly 0. It will be jittered before further use.")
    spca_conv_crit <- 1e-15
  }
  
  # Mishandling of check_rank
  parallel <- checkBoolean(parallel, "parallel")
  
  # Misshandling of svd_method
  if(!(svd_method %in% c("fast", "precise"))){
    stop(paste0("svd_method must be \"fast\" for usage of Eigen's BDCSVD or \"precise\" for usage of Eigen's JacobiSVD"))
  }
  
  # Misshandling of normalise
  normalise <- checkBoolean(normalise, "normalise")
  
  # Misshandling of comp_var_expl
  comp_var_expl <- checkBoolean(comp_var_expl, "comp_var_expl")
  
  result <- runSPCA(
    X_in = data_r,
    delay = delay, 
    selected = selected, 
    R = as.integer(no_of_factors), 
    ridge_penalty, 
    lasso_penalty, 
    max_iterations = max_iterations, 
    steps = max_no_steps, 
    weights = weights,
    comp_null = comp_null, 
    spca_conv_crit = spca_conv_crit, 
    parallel = parallel,
    svd_method = svd_method, 
    normalise = normalise, 
    comp_var_expl = comp_var_expl)
  
  # Rename the results
  names(result) <- c("loading_matrix_estim", "factor_estimate", "total_var_expl",
                     "pct_var_expl")
  result$data <- data
  
  if(is.zoo(data) || is.xts(data)){ # Also convert factors to time series
    start_vector <- c(year(time(data)[1]), month(time(data)[1]))
    result$factor_estimate <- as.zoo(ts(t(result$factor_estimate), start = start_vector, frequency = 12))
    colnames(result$factor_estimate) <- paste0("Factor ", 1:no_of_factors)
  }
  
  result$call <- func_call
  class(result) <- "SPCAFit"
  return(result)
}

#' @name print.SPCAFit
#' @title Generic printing function for SPCAFit S3 objects
#' @description
#' Print a compact summary of an `SPCAFit` object.
#'
#' @param x `SPCAFit` object.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' No return value; Prints a summary to the console.
#'
#' @author
#' Domenic Franjic
#' 
#' @export
print.SPCAFit <- function(x, ...) {
  simulated_time_series <- is.zoo(x$factor_estimate)
  no_of_factors <- ifelse(simulated_time_series, dim(x$factor_estimate)[2], dim(x$factor_estimate)[1])
  no_of_obs <- ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2])
  cat("Simulated Dynamic Factor Model\n")
  cat("=========================================================================\n")
  cat("No. of Observations                        :", ifelse(simulated_time_series, dim(x$data)[1], dim(x$data)[2]), "\n")
  cat("No. of Variables                           :", ifelse(simulated_time_series, dim(x$data)[2], dim(x$data)[1]), "\n")
  cat("No. of Factors                             :", no_of_factors, "\n")
  cat("No. of zero elements in the loading matrix :", sum(x$loading_matrix_estim == 0), "\n")
  cat("=========================================================================\n")
  cat("Head of the factors :\n")
  max_print <- min(5, no_of_obs)
  if(simulated_time_series){
    print(head(x$factor_estimate, max_print))
  }else{
    print(x$factor_estimate[, 1:max_print])
  }
  cat("Tail of the factors :\n")
  if(simulated_time_series){
    print(tail(x$factor_estimate, max_print))
  }else{
    print(x$factor_estimate[, (dim(x$factor_estimate)[2] - (max_print - 1)):(dim(x$factor_estimate)[2])])
  }
  max_print_loadings <- min(5, ifelse(simulated_time_series, dim(x$factor_estimate)[1], dim(x$factor_estimate)[2]))
  cat("Head of the loading matrix :\n")
  print(head(x$loading_matrix_estim, max_print_loadings))
  cat("Tail of the loading matrix :\n")
  print(tail(x$loading_matrix_estim, max_print_loadings))
  cat("=========================================================================\n")
  invisible(x)
}

#' @name plot.SPCAFit
#' @title Generic plotting function for SPCAFit S3 objects
#' @description
#' Create diagnostic plots for an `SPCAFit` object.
#'
#' @param x `SPCAFit` object.
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
#'   \item{`Factor Time Series Plots`}{`patchwork`/`ggplot` object showing the 
#'   estimated factors over time.}
#'   \item{`Loading Matrix Heatmap`}{`ggplot` object showing a heatmap of the 
#'   estimated factor loadings. Zeros are highlighted in black.}
#'   \item{`Meas. Error Var.-Cov. Matrix Heatmap`}{`ggplot` object showing a 
#'   heatmap of the measurement error variance–covariance matrix.}
#'   \item{`Eigenvalue Plot`}{`ggplot` object showing a bar plot of the 
#'   eigenvalues of the measurement error variance–covariance matrix.}
#'  }
#' 
#' @author
#' Domenic Franjic
#' 
#' @export
plot.SPCAFit <- function (x, 
                          axis_text_size = 20, 
                          legend_title_text_size = 20,
                          ...) {
  out_list <- list()
  if (is.zoo(x$data)) {
    series_names <- colnames(x$data)
    no_of_factors <- dim(x$factor_estimate)[2]
    no_of_obs <- dim(x$factor_estimate)[1]
    time_vector <- as.Date(time(x$factor_estimate))
    factors <- x$factor_estimate
  } else {
    series_names <- rownames(x$data)
    no_of_factors <- dim(x$factor_estimate)[1]
    no_of_obs <- dim(x$factor_estimate)[2]
    time_vector <- 1:dim(x$factor_estimate)[2]
    factors <- t(x$factor_estimate)
    factors <- as.zoo(ts(factors, start = c(1, 1), frequency = 12))
  }
  
  out_list$`Factor Time Series Plots` <- plotFactorEstimates(factors, matrix(0, no_of_factors, no_of_factors * no_of_obs), 
                                                             no_of_factors, axis_text_size)
  
  out_list$`Loading Matrix Heatmap` <- plotLoadingHeatMap(x$loading_matrix_estim, series_names, 
                                                          no_of_factors, axis_text_size, 
                                                          legend_title_text_size)
  
  if (is.zoo(x$data)) {
    residuals <- coredata(na.omit(x$data)) - coredata(x$factor_estimate) %*% 
      t(x$loading_matrix_estim)
  } else {
    residuals <- na.omit(t(x$data)) - t(x$factor_estimate) %*% 
      t(x$loading_matrix_estim)
  }
  measurement_error_var_cov_df <- as.data.frame(t(residuals) %*% residuals * 1/(dim(residuals)[1] - 1))
  out_list$`Meas. Error Var.-Cov. Matrix Heatmap` <- plotMeasVarCovHeatmap(measurement_error_var_cov_df, 
                                                                           series_names, axis_text_size, 
                                                                           legend_title_text_size)
  
  out_list$`Eigenvalue Plot` <- plotMeasVarCovEigenvalues(eigen(measurement_error_var_cov_df)$values, 
                                                          no_of_factors, axis_text_size, 
                                                          legend_title_text_size)
  
  if (length(x$pct_var_expl) != 0) {
    var_explained_df <- data.frame(Component = c(paste0("Factor ", 1:no_of_factors), "Unexplained"), 
                                   `Variance Explained` = c(x$pct_var_expl, 1 - sum(x$pct_var_expl)), 
                                   check.names = FALSE)
    var_explained_df$Component <- factor(var_explained_df$Component)
    colourPalette <- grDevices::colorRamp(c("#88ccee", "#FFFFFF", "#117733"))
    cols_mat <- colourPalette(seq(0, 1, length.out = no_of_factors + 1))
    pie_chart_colours <- grDevices::rgb(cols_mat[, 1], cols_mat[, 2], cols_mat[, 3], maxColorValue = 255)
    names(pie_chart_colours) <- levels(var_explained_df$Component)
    var_explained_df$pct_label <- sprintf("%.3f%%", var_explained_df$`Variance Explained`)
    max_radial_position <- 1 - cumsum(var_explained_df$`Variance Explained`)
    min_radial_position <- c(1, head(max_radial_position, -1))
    var_explained_df$radial_position <- (min_radial_position + max_radial_position)/2
    out_list$`Variance Explained Chart` <- ggplot(var_explained_df, aes(x = "", y = `Variance Explained`, fill = Component)) +
      geom_col(width = 1) + coord_polar(theta = "y") +
      scale_fill_manual(values = pie_chart_colours) + 
      geom_text(aes(y = radial_position, label = pct_label), color = "black", size = axis_text_size * 0.5) +
      theme_void() + theme(text = element_text(size = axis_text_size),
                           legend.title = element_text(size = legend_title_text_size)
                           )
  } else {
    out_list$`Variance Explained Chart` <- "The relative variance explained by each factor has not been computed as comp_var_expl = FALSE"
  }
  return(out_list)
}