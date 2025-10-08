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

install.packages("../TwoStepSDFM_0.0.0.2.tar.gz", repos = NULL, type = "source")

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
no_of_observations <- 200 # Number of observations
no_of_variables <- 20 # Number of variabes
no_of_factors <- 1 # Number of factors
trans_error_var_cov <- diag(1, no_of_factors) # Variance-covariance matrix of the transition errors
loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors) # Factor loadings matrix
meas_error_mean <- rep(0, no_of_variables) # Mean of the measurement error
meas_error_var_cov <- diag(1, no_of_variables) # Variance-covariance matrix of the measurement error
trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors)) # Factor VAR coefficient matrix
factor_lag_order <- 2 # Order of the factor VAR process
simul_delay <- c(6, floor(rexp(19, 1)))
quarterfy <- TRUE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
quarterly_variable_ratio  <- 1/20
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
FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
            no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
            meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
            trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
            factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
            quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
            beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
            rescale = rescale, check_stationarity = check_stationarity, 
            stationarity_check_threshold = stationarity_check_threshold)

print(FM)
plots <- plot(FM)
plots$`Factor Time Series Plots`
plots$`Loading Matrix Heatmap`
plots$`Meas. Error Var.-Cov. Matrix Heatmap`

frequency <- FM$frequency
data <- FM$data[, which(frequency == 12)]
delay <- simul_delay[which(frequency == 12)]
selected <- c(sum(frequency == 12))
no_of_factors <- 1
max_factor_lag_order  <- 10
decorr_errors <- TRUE 
lag_estim_criterion  <- "BIC"
ridge_penalty <- 1e-06
lasso_penalty <- NULL
max_iterations <- 1000
max_no_steps <- NULL
comp_null <- 1e-15
check_rank <- FALSE
conv_crit <- 1e-04 
conv_threshold <- 1e-04
log <- FALSE
parallel <- FALSE
fcast_horizon <- 0

fit_sparse <- twoStepSDFM(data = data, delay = delay, selected = selected, 
                          no_of_factors = no_of_factors,  
                          max_factor_lag_order  = max_factor_lag_order, 
                          decorr_errors = decorr_errors, 
                          lag_estim_criterion  = lag_estim_criterion,
                          ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                          max_iterations = max_iterations, max_no_steps = max_no_steps, 
                          comp_null = comp_null, check_rank = check_rank,  
                          conv_crit = conv_crit, conv_threshold = conv_threshold, 
                          log = log, parallel = parallel, fcast_horizon = fcast_horizon)

fit_sparse$loading_matrix_estimate
fit_sparse$smoothed_factors
fit_sparse$smoothed_state_variance
fit_sparse$factor_var_lag_order
fit_sparse$error_var_cov_cholesky_factor
fit_sparse
print(fit_sparse)
plots <- plot(fit_sparse)
plots$`Factor Time Series Plots`
plots$`Loading Matrix Heatmap`
plots$`Meas. Error Var.-Cov. Matrix Heatmap`

data <- FM$data
variables_of_interest <- 1
max_fcast_horizon <- 4
delay <- simul_delay
frequency <- FM$frequency
selected <- c(sum(frequency == 12))
no_of_factors <- 1
max_factor_lag_order  <- 10
decorr_errors <- TRUE 
lag_estim_criterion  <- "BIC"
ridge_penalty <- 1e-06
lasso_penalty <- NULL
max_iterations <- 1000
max_no_steps <- NULL
comp_null <- 1e-15
check_rank <- FALSE
conv_crit <- 1e-04 
conv_threshold <- 1e-04
log <- FALSE
parallel <- FALSE
fcast_horizon <- 10
max_ar_lag_order <- 5
max_predictor_lag_order  <- 5

results <- nowcast(data = data, variables_of_interest = variables_of_interest, 
                   max_fcast_horizon = max_fcast_horizon, delay = delay, selected = selected,
                   frequency = frequency, no_of_factors = no_of_factors, 
                   max_factor_lag_order  = max_factor_lag_order,  
                   decorr_errors = decorr_errors, 
                   lag_estim_criterion  = lag_estim_criterion,
                   ridge_penalty = ridge_penalty, lasso_penalty = lasso_penalty, 
                   max_iterations = max_iterations, max_no_steps = max_no_steps, 
                   comp_null = comp_null, check_rank = check_rank,  
                   conv_crit = conv_crit, conv_threshold = conv_threshold, 
                   log = log, parallel = parallel, max_ar_lag_order = max_ar_lag_order,
                   max_predictor_lag_order = max_predictor_lag_order
)

results$Forecasts
results$`Single Predictor Forecasts`
results$`SDFM Fit`


print(results)
plots <- plot(results)
plots$`Single Pred. Fcast Density Plots Series 1`
plots$`Single Pred. Fcast Density Plots Series 2`
plots$`Single Pred. Fcast Density Plots Series 3`
plots$`Single Pred. Fcast Density Plots Series 4`
plots$`Single Pred. Fcast Density Plots Series 5`
plots$`Factor Time Series Plots`
plots$`Loading Matrix Heatmap`
