#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import zoo
#' @import xts
#' @import lubridate
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

#' @name noOfFactors
#' @title Estimate the number of Factors
#' @description
#' Estimate the number of factors of a linear Gaussian latent factor model using 
#' the eigenvalue slope test of \insertRef{onatski2009testing}{TwoStepSDFM} and 
#' the infromation criterion based approach of 
#' \insertRef{bai2002determining}{TwoStepSDFM}.
#' @param data Numeric (no_of_vars \eqn{\times}{x} no_of_obs) matrix of data or 
#' zoo/xts object sampled at the same frequency.
#' @param min_no_factors Integer minimum number of factors to be tested.
#' @param max_no_factors Integer maximum number of factors to be tested (should 
#' be at most min_no_factors + 17).
#' @param confidence_threshold Numeric threshold value to stop the testing procedure.
#' 
#' @details
#' The \insertCite{onatski2009testing}{TwoStepSDFM} procedure splits the data 
#' matrix along the time dimension into two equally sized (`no_of_vars` 
#' \eqn{\times}{x} `cut_off`) sub-matrices \eqn{\bm{X}_{1/2}}{`data_fst_half`} 
#' and \eqn{\bm{X}_{2/2}}{`data_snd_half`}. It then proceeds to build 
#' \eqn{\tilde{\bm{X}} := \bm{X}_{1/2} + i\bm{X}_{2/2}}{`complex_data = data_fst_half + i * data_snd_half`},
#' where \eqn{i=\sqrt{-1}}{`i=sqrt(-1)`}. We then compute the eigenvalues of the 
#' Gram matrix \eqn{\tilde{\bm{X}} 
#' \tilde{\bm{X}}^{\dagger}}{`complex_data %*% Conj(t(complex_data))`}, where 
#' \eqn{\tilde{\bm{X}}^{\dagger}}{`Conj(t(complex_data))`} represents the 
#' adjoint. Finally, a test based on the computed eigenvalues is performed. 
#' This test is an iterative testing procedure, starting by testing the null
#' that the true number of factors is `min_no_factors`. If the test is rejected
#' by comparison of the \eqn{p}{p}-value against `confidence_threshold`, we 
#' test whether the true number of factors is `min_no_factors + 1` until we can 
#' no longer reject at `confidence_threshold` or `max_no_factors` is reached.
#' 
#' As the distribution of the eigenvalues under the null is nonstandard
#' \insertCite{onatski2009testing}{TwoStepSDFM}, 
#' simulated critical values are used. They are retrieved from
#' \insertRef{onatski2009testing_supl}{TwoStepSDFM}. As the range of the 
#' simulated critical values is limited, the minimum and maximum number of 
#' potential factors is limited such that `max_no_factors` should be no more 
#' than `min_no_factors + 17`. However, it is recommended to operate well below 
#' this maximum as the test size decreases with 
#' `max_no_factors - min_no_factors`. 
#' 
#' The \insertCite{bai2002determining}{TwoStepSDFM} information criterion 
#' determines the number of factors by minimising a BIC. Here, three different 
#' penalty terms are provided. It is up to the user to determine the most 
#' appropriate for the problem at hand. In general, however, the second 
#' information criterion is used.
#' 
#' @return 
#' An object of class `NoOfFactorsFit` with components:
#' \describe{
#'   \item{no_of_factors}{Integer estimated number of factors.}
#'   \item{p_value}{Numeric \eqn{p}{p}-value of the final test.}
#'   \item{confidence_threshold}{Numeric significance level used.}
#'   \item{statistic}{Numeric test statistic value of the last test.}
#'   \item{eigen_values}{Numeric vector of eigenvectors of the complex data Gram
#'   matrix.}
#' }
#'
#' @author
#' Domenic Franjic
#' 
#' @references
#' \insertRef{onatski2009testing}{TwoStepSDFM}
#' 
#' \insertRef{onatski2009testing_supl}{TwoStepSDFM}
#'
#' @examples
#' data(factor_model)
#' no_of_factors_estim <- noOfFactors(data = factor_model$data, min_no_factors = 1, 
#'                                    max_no_factors = 5, confidence_threshold = 0.05)
#' print(no_of_factors_estim)
#' factor_estim_plots <- plot(no_of_factors_estim)
#' factor_estim_plots$`Eigen Value Plot Test Procedure`
#' factor_estim_plots$`IC plot for IC1`
#' factor_estim_plots$`IC plot for IC2`
#' factor_estim_plots$`IC plot for IC3`
#' 
#' @export
noOfFactors <- function(data, min_no_factors = 1, max_no_factors = 7, confidence_threshold = 0.05){
  func_call <- match.call()
  
  # Mishandling of data
  if(!is.zoo(data) && !is.xts(data)){
    data_r <- try(as.matrix(data), silent = TRUE)
    if (inherits(data_r, "try-error")) {
      stop(paste0("data must be a matrix, convertible to a matrix or a time-series/zoo object"))
    }
  }else{
    data_r <- try(t(coredata(data)), silent = TRUE)
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
  na_ind <- -unique(which(is.na(data_r), arr.ind = TRUE)[, 2])
  if(length(na_ind) != 0){
    message(paste0("Cut ", length(na_ind)," observations due to NAs."))
    no_na_data <- as.matrix(data_r[, na_ind, drop = FALSE])
  }else{
    no_na_data <- as.matrix(data_r[, , drop = FALSE])
  }
  
  
  # Mishandling of max_no_factors and min_no_factors
  max_no_factors <- checkPositiveSignedInteger(max_no_factors, "max_no_factors");
  if(max_no_factors >= dim(no_na_data)[1] - 2){
    stop(paste0("max_no_factors must be smaller than dim(no_na_data)[1] - 2 = ", dim(no_na_data)[1] - 2, "."))
  }
  min_no_factors <- checkPositiveSignedInteger(min_no_factors, "min_no_factors");
  if(min_no_factors <= 0){
    stop(paste0("min_no_factors must be strictly positive."))
  }
  if(min_no_factors >= max_no_factors){
    stop(paste0("max_no_factors must be strictly greater than min_no_factors."))
  }
  if (7 < max_no_factors - min_no_factors) {
    warning(paste0("Power of the test might be low as max_no_factors - min_no_factors = ", 
                   max_no_factors - min_no_factors," > 7."))
  }
  if (18 < max_no_factors - min_no_factors) {
    stop(paste0("Critical values for max_no_factors - min_no_factors = ", max_no_factors - min_no_factors, 
                " > 18 not available. Decrease max_no_factors"))
  }
  
  # Mishandling of confidence_threshold = 0.05
  confidence_threshold <- checkPositiveDouble(confidence_threshold, "confidence_threshold")
  if(confidence_threshold <= 0 || confidence_threshold >= 1){
    stop(paste0("confidence_threshold must be in (0,1)."))
  }
  
  # The values for the test-statistics stem: https://www.econometricsociety.org/publications/econometrica/2009/09/01/testing-hypotheses-about-number-factors-large-factor-models (Last accessed: 25.11.2025, 10:03)
  file_path <- system.file("extdata", "Onatski_test_stats_csv.txt", package = "TwoStepSDFM")
  test_values <- as.matrix(read.table(file_path, sep = ",", header = FALSE))
  results <- list()
  results$test <- runNoOfFactorsTest(no_na_data, test_values, min_no_factors, max_no_factors, confidence_threshold)
  results$ic <- runNoOfFactorsInfoCrit(no_na_data, max_no_factors)
  results$ic$ic_one_min <- which.min(results$ic$information_crit[, 1])
  results$ic$ic_two_min <- which.min(results$ic$information_crit[, 2])
  results$ic$ic_three_min <- which.min(results$ic$information_crit[, 3])
  results$max_no_factors <- max_no_factors
  
  results$call <- call
  class(results) <- "NoOfFactorsFit"
  return(results)
}


#' @name print.NoOfFactorsFit
#' @title Generic printing function for NoOfFactorsFit S3 objects
#' @description
#' Print a compact summary of an `NoOfFactorsFit` object.
#'
#' @param x `NoOfFactorsFit` object.
#' @param ... Additional parameters for the plotting functions.
#'
#' @return
#' No return value; Prints a summary to the console.
#'
#' @author
#' Domenic Franjic
#' 
#' @export
print.NoOfFactorsFit <- function(x, ...) {
  cat("Results of the Onatsky testing procedure:\n")
  cat(paste0("The estimated no. of factors is ", x$test$no_of_factors, " with a p-value of ", x$test$p_value, " and a critical value of alpha = ", x$test$confidence_threshold))
  if(x$test$no_of_factors == x$max_no_factors - 1){
    cat("\n")
    warning(paste0("No. of factors has been chosen as max_no_factors - 1 = ", x$max_no_factors - 1, 
                   ". It might be necessary to increase max_no_factors and repeat the procedure"))
  }
  cat("\n\n")
  cat("Results of the Bai and Ng information criteria:\n")
  cat(paste0("No. of factors according to IC1: ", x$ic$ic_one_min, "\n"))
  cat(paste0("No. of factors according to IC2: ", x$ic$ic_one_min, "\n"))
  cat(paste0("No. of factors according to IC3: ", x$ic$ic_one_min, "\n"))
  invisible(x)
}

#' @name plot.NoOfFactorsFit
#' @title Generic plotting function for NoOfFactorsFit S3 objects
#' @description
#' Create diagnostic plots for an `NoOfFactorsFit` object, 
#'
#' @param x `NoOfFactorsFit` object.
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
#'   \item{`Eigen Value Plot`}{`ggplot` object showing a bar plot of the 
#'   eigenvalues of the complex data Gram matrix.}
#' }
#'
#' @author
#' Domenic Franjic
#' 
#' @export
plot.NoOfFactorsFit <- function(x, 
                                axis_text_size = 20, 
                                legend_title_text_size = 20, 
                                ...) {
  out_list <- list()
  
  # Result plot for Onatsky testing procedure
  out_list$`Eigen Value Plot Test Procedure` <- plotMeasVarCovEigenvalues(x$test$eigen_values,
                                                                          x$test$no_of_factors,
                                                                          axis_text_size,
                                                                          legend_title_text_size)
  
  # Result plots for IC criteria of Bai and Ng
  out_list$`IC plot for IC1` <- plotInformationCrit(x$ic$information_crit[, 1], x$ic$ic_one_min, axis_text_size)
  out_list$`IC plot for IC2` <- plotInformationCrit(x$ic$information_crit[, 2], x$ic$ic_two_min, axis_text_size)
  out_list$`IC plot for IC3` <- plotInformationCrit(x$ic$information_crit[, 3], x$ic$ic_three_min, axis_text_size)
  
  return(out_list)
}




