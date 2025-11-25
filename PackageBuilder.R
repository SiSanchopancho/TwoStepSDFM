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
# devtools::install()
# devtools::check()
roxygen2::roxygenise()

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

output <- system2("R", args = c("CMD", "INSTALL", "../TwoStepSDFM_0.1.3.tar.gz"),
                  stdout = TRUE, stderr = TRUE)
writeLines(output, "compilation_log.txt")
tail(output)

# install.packages("../TwoStepSDFM_0.1.1.tar.gz", repos = NULL, type = "source")

library(TwoStepSDFM)
library(testthat)
library(zoo)
library(xts)
library(lubridate)
library(ggplot2)
ls("package:TwoStepSDFM")

set.seed(02102025)
no_of_observations <- 200
no_of_variables <- 10
no_of_factors <- 3
trans_error_var_cov <- diag(1, no_of_factors)
loading_matrix <- matrix(round(rnorm(no_of_variables * no_of_factors)), no_of_variables, no_of_factors)
meas_error_mean <- rep(0, no_of_variables)
meas_error_var_cov <- diag(1, no_of_variables)
trans_var_coeff <- cbind(diag(0.5, no_of_factors), -diag(0.25, no_of_factors))
factor_lag_order <- 2 
simul_delay <- c(floor(rexp(10, 1)))
quarterfy <- FALSE
quarterly_variable_ratio  <- 0
corr <- TRUE
beta_param <- 2
seed <- 01102025
set.seed(seed)
burn_in <- 999
starting_date <- "1970-01-01"
rescale <- TRUE
check_stationarity <- TRUE
stationarity_check_threshold <- 1e-10

FM <- simFM(no_of_observations = no_of_observations, no_of_variables = no_of_variables, 
            no_of_factors = no_of_factors, loading_matrix = loading_matrix, 
            meas_error_mean = meas_error_mean, meas_error_var_cov = meas_error_var_cov,
            trans_error_var_cov = trans_error_var_cov, trans_var_coeff = trans_var_coeff, 
            factor_lag_order = factor_lag_order, delay = simul_delay, quarterfy = quarterfy, 
            quarterly_variable_ratio  = quarterly_variable_ratio, corr = corr, 
            beta_param = beta_param, seed = seed, burn_in = burn_in, starting_date = starting_date,
            rescale = rescale, check_stationarity = check_stationarity, 
            stationarity_check_threshold = stationarity_check_threshold)

frequency <- FM$frequency
delay <- simul_delay
no_of_factors <- 3 
max_factor_lag_order  <- 10
decorr_errors <- TRUE 
lag_estim_criterion  <- "BIC"
comp_null <- 1e-15
check_rank <- FALSE
log <- FALSE
parallel <- FALSE
fcast_horizon <- 10

fit <- twoStepDenseDFM(data = FM$data, delay = delay, no_of_factors = no_of_factors,
                       max_factor_lag_order  = max_factor_lag_order, 
                       decorr_errors = decorr_errors,  lag_estim_criterion  = lag_estim_criterion, 
                       comp_null = comp_null,  check_rank = check_rank, log = log, 
                       parallel = parallel, fcast_horizon = fcast_horizon)

print(fit)
graphs <- plot(fit)
graphs$`Factor Time Series Plots`
graphs$`Loading Matrix Heatmap`
graphs$`Meas. Error Var.-Cov. Matrix Heatmap`