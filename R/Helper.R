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

#' Helper function to check unsigned integer function parameter
#' @keywords internal
checkPositiveSignedInteger <- function(parameter, name, bit_size = 64) {
  if(!is.numeric(parameter)){
    stop(paste0(name, " is not numeric."))
  }
  if(length(parameter) != 1 || is.na(parameter)) {
    stop(paste0(name, " must be a single, non-NA numeric value."))
  }
  if(is.infinite(parameter)){
    stop(paste0(name, " cannot be Inf."))
  }
  if(parameter %% 1 != 0){
    warning(paste0(name, " is not an integer. It will be truncated before further use."))
    parameter <- floor(parameter)
  }
  if(parameter < 0 || parameter > 2^(bit_size - 1) - 1) {
    stop(paste0(name, " must be a non-negative signed 64bit integer."))
  }
  
  return(parameter)
}

#' Helper function to check parameter matrices
#' @keywords internal
checkParameterMatrix <- function(parameter, name, no_of_rows, no_of_cols){
  parameter <- try(as.matrix(parameter), silent = TRUE)
  if(inherits(parameter, "try-error")) {
    stop(paste0(name, " must be a matrix or convertible to a matrix."))
  }
  if(!is.numeric(parameter)){
    stop(paste0(name, " has non-numeric elements."))
  }
  if (dim(parameter)[1] != no_of_rows || dim(parameter)[2] != no_of_cols) {
    stop(paste0(name, " must be of dimensions (", no_of_rows, "x", no_of_cols, ") but is (", dim(parameter)[1], "x", dim(parameter)[2], ")."))
  }
  if(any(is.infinite(parameter))){
    stop(paste0(name, " cannot have (-)Inf values."))
  }
  if(any(is.na(parameter))){
    stop(paste0(name, " has NA values."))
  }

  return(parameter)
}

#' Helper function to check parameter vectors for positive signed integer values
#' @keywords internal
checkPositiveSignedParameterVector <- function(parameter, name, size){
  parameter <- try(as.matrix(parameter), silent = TRUE)
  if(inherits(parameter, "try-error")) {
    stop(paste0(name, " must be a matrix or convertible to a matrix"))
  }
  if(all(dim(parameter) == c(1, size))) {
    parameter <- t(parameter)
  }
  if(any(dim(parameter) != c(size, 1))) {
    stop(paste0(name, " must be a vector of a matrix of dimensions (", size, 
                "x1) or (1x", size, "). The provided matrix is (", dim(parameter)[1], 
                "x", dim(parameter)[2], ")."))
  }
  if(!is.numeric(parameter)){
    stop(paste0(name, " has non-numeric elements."))
  }
  if(any(is.na(parameter))){
    stop(paste0(name, " has NA values."))
  }
  if(any(parameter %% 1 != 0)) {
    warning(paste0("At least one element of ", name, " is not an integer. It will be truncated before further use."))
    parameter <- floor(parameter)
  }
  if(any(parameter < 0) || any(parameter > 2^63 - 1)) {
    stop(paste0(name, " must be a non-negative signed 64bit integer."))
  }

  return(parameter)
}

#' Helper function to check parameter vectors for positive signed integer values
#' @keywords internal
checkBoolean <- function(parameter, name){
  if(is.null(parameter))
  {
    stop(paste0(name, " must be a single, non-NA boolean value."))
  }
  if(!is.logical(parameter) && !(parameter %in% c(0, 1))){
    stop(paste0(name, " must be a single, non-NA boolean value."))
  }
  if (length(parameter) != 1 || is.na(parameter)) {
    stop(paste0(name, " must be a single, non-NA boolean value."))
  }
  
  return(as.logical(parameter))
}

#' Helper function to check positive double function parameter
#' @keywords internal
checkPositiveDouble <- function(parameter, name) {
  if(!is.numeric(parameter)){
    stop(paste0(name, " is not numeric."))
  }
  if(is.infinite(parameter)){
    stop(paste0(name, " cannot be (-)Inf."))
  }
  if (length(parameter) != 1 || is.na(parameter)) {
    stop(paste0(name, " must be a single, non-NA numeric value."))
  }
  if (parameter < 0) {
    stop(paste0(name, " must be non-negative."))
  }
  
  return(parameter)
}

#' Helper function to check positive double function parameter
#' @keywords internal
makeRaggedEdges <- function(data, delay){
no_of_observations <- dim(data)[1]
no_of_variables <- dim(data)[2]
for(n in 1:no_of_variables){
  if(delay[n] == 0){
    next
  }
  data[(no_of_observations + 1 - delay[n]):(no_of_observations), n] <- NaN
}
return(data)
}


