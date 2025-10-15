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

# output <- system2("R", args = c("CMD", "INSTALL", "../TwoStepSDFM_0.0.0.2.tar.gz"), 
# stdout = TRUE, stderr = TRUE)
# writeLines(output, "compilation_log.txt")

install.packages("../TwoStepSDFM_0.1.1.tar.gz", repos = NULL, type = "source")

library(TwoStepSDFM)
library(testthat)
library(zoo)
library(xts)
library(lubridate)
library(ggplot2)
ls("package:TwoStepSDFM")

# Test

# Simulate a DGP using simFM
set.seed(02102025)
no_of_observations <- 300 + 3# Number of observations
no_of_variables <- 50 # Number of variabes
no_of_factors <- 2 # Number of factors
trans_error_var_cov <- diag(1, no_of_factors) # Variance-covariance matrix of the transition errors
loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors) # Factor loadings matrix
meas_error_mean <- rep(0, no_of_variables) # Mean of the measurement error
meas_error_var_cov <- diag(1, no_of_variables) # Variance-covariance matrix of the measurement error
trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors)) # Factor VAR coefficient matrix
factor_lag_order <- 2 # Order of the factor VAR process
simul_delay <- c(3, floor(rexp(no_of_variables - 1, 1)))
quarterfy <- TRUE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
quarterly_variable_ratio  <- 1/no_of_variables
corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
beta_param <- 2 # Beta parameter governing the degree of correlation of the measurement error
seed <- 01102025 # Seed
set.seed(seed)
burn_in <- 999 # Burn-in period
starting_date <- "1970-01-01"
rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 
check_stationarity <- TRUE
stationarity_check_threshold <- 1e-10

# Draw the FM object
is_FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
               no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
               meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
               trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
               factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
               quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
               beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
               rescale = rescale, check_stationarity = check_stationarity, 
               stationarity_check_threshold = stationarity_check_threshold)

# Cross validation loop
library(doSNOW)
library(doParallel)
library(foreach)

data <- is_FM$data
variable_of_interest <- 1
fcast_horizon <- 0
delay <- simul_delay
frequency <- is_FM$frequency
no_of_factors <- 3
seed <- 09102025
min_ridge_penalty <- 0.01
max_ridge_penalty <- 1
lasso_penalty_type <- "selected"
min_max_penalty <- c(10, no_of_variables - 1)
# min_max_penalty <- c(10, 50 * (no_of_variables - 1))
# min_max_penalty <- c(0.0001, 10)
cv_repititions <- 3
cv_size <- 1000
max_factor_lag_order = 10
decorr_errors = TRUE
lag_estim_criterion = "BIC"
ridge_penalty = 1e-6
lasso_penalty = NULL
max_iterations = 5000
max_no_steps = NULL
comp_null = 1e-15
check_rank = FALSE
conv_crit = 1e-4
conv_threshold = 1e-4
log = FALSE
parallel = TRUE
max_ar_lag_order = 5
max_predictor_lag_order = 5

cv <- crossVal(data = data, variable_of_interest = variable_of_interest, fcast_horizon = fcast_horizon,
               delay = delay, frequency = frequency, no_of_factors = no_of_factors,
               seed = seed, min_ridge_penalty = min_ridge_penalty, max_ridge_penalty = max_ridge_penalty,
               cv_repititions = cv_repititions, cv_size = cv_size, lasso_penalty_type = lasso_penalty_type,
               min_max_penalty = min_max_penalty, max_factor_lag_order = max_factor_lag_order,
               decorr_errors = decorr_errors, lag_estim_criterion = lag_estim_criterion,
               ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, max_iterations = max_iterations,
               max_no_steps = max_no_steps, comp_null = comp_null, check_rank = check_rank,
               conv_crit = conv_crit, conv_threshold = conv_threshold, parallel = parallel,
               max_ar_lag_order = max_ar_lag_order, max_predictor_lag_order = max_predictor_lag_order)
cv
a <- plot(cv)
a$`CV Results`
a$`BIC Results`







