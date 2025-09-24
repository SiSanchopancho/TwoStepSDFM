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

# Convinient package builder #

# Load libraries
library(rstudioapi)
library(roxygen2)

# Set directory
setwd(dirname(getActiveDocumentContext()$path))
rm(list = ls())

# Clean up
unlink("src/*.o")
unlink("src/*.so")
unlink("src/*.dll")
devtools::clean_dll()
devtools::clean_vignettes()
devtools::build_vignettes()

# Compile Rcpp attributes
Rcpp::compileAttributes()

# Install the package
#devtools::install()
devtools::check()

# Build
devtools::build()

detach("package:TwoStepSDFM", unload = TRUE)

.rs.restartR()

# Load libraries
library(rstudioapi)
library(roxygen2)

# Set directory
setwd(dirname(getActiveDocumentContext()$path))
rm(list  = ls())

install.packages("../TwoStepSDFM_0.0.0.2.tar.gz", repos = NULL, type = "source")

library(TwoStepSDFM)
library(testthat)
ls("package:TwoStepSDFM")
# Test

# Simulate a DGP using simFM
no_of_observations <- 100 # Number of observations
no_of_variables <- 10 # Number of variabes
no_of_factors <- 3 # Number of factors
trans_error_var_cov <- diag(1, no_of_factors) # Variance-covariance matrix of the transition errors
loading_matrix <- matrix(rnorm(no_of_variables * no_of_factors), no_of_variables, no_of_factors) # Factor loadings matrix
meas_error_mean <- rep(0, no_of_variables) # Mean of the measurement error
meas_error_var_cov <- diag(1, no_of_variables) # Variance-covariance matrix of the measurement error
trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors), -diag(0.25, no_of_factors), -diag(0.25, no_of_factors), -diag(0.25, no_of_factors)) # Factor VAR coefficient matrix
factor_lag_order <- 5 # Order of the factor VAR process
quarterfy <- FALSE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
quarterly_variable_ratio  <- 0
corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
seed <- 16022024 # Seed
set.seed(seed)
burn_in <- 999 # Burn-in period
rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 

# Draw the FM object
FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, no_of_factors = no_of_factors, loading_matrix = loading_matrix, meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
            trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, factor_lag_order = factor_lag_order, quarterfy = quarterfy, quarterly_variable_ratio  = quarterly_variable_ratio ,
            corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale,
            check_stationarity = TRUE, stationarity_check_threshold = 1e-10)

# Fitting a sparse model with l2 regularisation and non-orthogonal measurement errors
selected <- c(round(no_of_variables * 0.8), round(no_of_variables * 0.5), round(no_of_variables * 0.5))
delay <- rep(0, no_of_variables)
fit_sparse <- twoStepSDFM(data = FM$data, delay = delay, selected = selected, no_of_factors = no_of_factors, 
                          max_factor_lag_order  = 10, decorr_errors = TRUE, 
                          lag_estim_criterion  = "BIC", ridge_penalty = 1e-06, 
                          lasso_penalty = NaN, max_iterations = 1000, max_no_steps = NaN, 
                          comp_null = 1e-15,  check_rank = FALSE,  conv_crit = 1e-04, 
                          conv_threshold = 1e-04, log = FALSE, parallel = FALSE)
fit_sparse$factor_estimate


