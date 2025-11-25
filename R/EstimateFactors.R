#' @useDynLib TwoStepSDFM, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import zoo
#' @import xts
#' @import lubridate
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

#' @name noOfFactors
#' @title Estimate the number of factors using the Onatski (2010) procedure.
#' @param data Numeric (no_of_variables x no_of_observations) matrix of data or zoo/xts object.
#' @param min_no_factors Integer minimum number of factors to be tested.
#' @param max_no_factors Integer maximum number of factors to be tested (should be at most min_no_factors + 17).
#' @param confidence_threshold Numeric threshold value to stop the testing procedure.
#' @return Returns a list with the number of factors and some additional information.
#' @export
noOfFactors <- function(data, min_no_factors = 1, max_no_factors = 15, confidence_threshold = 0.05){
  
  # Mishandling of data
  # Misshandling of the data matrix
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
    print(paste0("Cut ", length(na_ind)," observations due to NAs."))
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
  results <- runNoOfFactors(no_na_data, test_values, min_no_factors, max_no_factors, confidence_threshold)
  
  if(results$no_of_factors == max_no_factors - 1){
    warning(paste0("No. of factors has been chosen as max_no_factors - 1 = ", max_no_factors - 1, 
                   ". It might be necessary to increase max_no_factors and repeat the procedure"))
  }
  
  print(paste0("The estimated no. of factors is ", results$no_of_factors, " with a p-value of ", results$p_value, " and a critical value of alpha = ", results$confidence_threshold))
  
  return(results)
}
