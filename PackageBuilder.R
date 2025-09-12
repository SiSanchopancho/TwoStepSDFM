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

# Clean up
unlink("src/*.o")
unlink("src/*.so")
unlink("src/*.dll")
devtools::clean_dll()
devtools::clean_vignettes()

# Compile Rcpp attributes
Rcpp::compileAttributes()

# Create Rd + NAMESPACE
devtools::document()

# Install the package
devtools::build()

detach("package:TwoStepSDFM", unload = TRUE)

.rs.restartR()

install.packages("../TwoStepSDFM_0.0.0.2.tar.gz", repos = NULL, type = "source")

library(TwoStepSDFM)
ls("package:TwoStepSDFM")
# Test

# Simulate a DGP using simFM
T <- 600 # Number of observations
N <- 10 # Number of variabes
R <- 2 # Number of factors
Sigma_epsilon <- diag(1, R) # Variance-covariance matrix of the transition errors
Lambda <- matrix(rnorm(N * R), N, R) # Factor loadings matrix
mu_xi <- rep(0, N) # Mean of the measurement error
Sigma_xi <- diag(1, N) # Variance-covariance matrix of the measurement error
Phi <- cbind(diag(0.5, R), -diag(0.25, R)) # Factor VAR coefficient matrix
P <- 2 # Order of the factor VAR process
quarterfy <- FALSE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
m <- 0
corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
seed <- 16022024 # Seed
set.seed(seed)
burn_in <- 999 # Burn-in period
rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 

# Draw the FM object
FM <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
            Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = quarterfy, m = m,
            corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale,
            check_staionarity = TRUE, stationarity_check_threshold = 1e-10)

# Fitting a sparse model with l2 regularisation and non-orthogonal measurement errors
selected <- c(round(N * 0.8), round(N * 0.5))
delay <- round(runif(N, 0, 10))
fit_sparse <- twoStepSDFM(FM$X, delay, selected, R, l2 = 1e+4)


