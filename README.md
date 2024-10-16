# TwoStepSDFM
A ``C++``-based ``R`` implementation of a two-step estimation procedure for a (linear Gaussian) Sparse Dynamic Factor Model (SDFM) as outlined in Franjic and Schweikert (2024).

## Introduction

The ``TwoStepSDFM`` package provides a fast implementation of the Kalman Filter and Smoother (hereinafter KFS, see Koopman and Durbin, 2000) to estimate factors in a mixed-frequency SDFM framework, explicitly accounting for cross-sectional correlation in the measurement error. The KFS is initialized using results from Sparse Principal Components Analysis (SPCA) by Zou and Hastie (2006) in a preliminary step. This approach generalizes the two-step estimator for approximate dynamic factor models by Giannone, Reichlin, and Small (2008) and Doz, Giannone, and Reichlin (2011).
We estimate the loadings matrix Λ and factors fₜ of the following state-space model:

xₜ = Λ fₜ + ξₜ

fₜ = Σₚ Φₚ fₜ₋ₚ + εₜ

for t = 1, ..., T. We use SPCA to estimate Λ, resulting in exact zero elements due to regularization. We also explicitly model the non-diagonal variance-covariance matrix of the measurement error ξₜ (for more details see Franjic and Schweikert (2024)).

## Features

As this is an early beta version of the package, only a limited set of functions is currently available. Additional functions, such as validation wrappers, may be added in future updates.

- **Fast Model Simulation**: The ``simFM()`` function provides a flexible framework to simulate approximate DFMs.
- **Fast Model Estimation**: The ``twoStepSDFM()`` function provides a fast and convenient implementation of the two-step estimator outlined in Franjic and Schweikert (2024).
- **Compatibility**: All functions take advantage of ``C++`` for enhanced speed.
- **Open-Source**: Distributed under the GNU General Public License v3.0.

## Prerequisites

- **Eigen** (version 3.4.0 or later): A `C++` template library for linear algebra [@eigenweb]. [Eigen Website](https://eigen.tuxfamily.org/)
- **Rcpp**: A package for integrating `C++` code into `R` [@Eddelbuettel2011Rcpp]. [Rcpp CRAN repository](https://cran.r-project.org/web/packages/Rcpp/index.html)
- **RcppEigen**: A package for integrating the `Eigen` linear algebra library into `R` [@Bates2013EcppEigen]. [RcppEigen CRAN repository](https://cran.r-project.org/web/packages/RcppEigen/index.html)

## Installation

Currently, the package source files come with a zipped version of ``Eigen3``. `Rcpp` and `RcppEigen` can be downloaded from CRAN or directly installed from within `R`by calling ``install.packages("...")``.

To install the package itself, a short `R` script is provided (see `PackageBuilder.R`). The routine simply calls 
```{R, eval=FALSE}
system("R CMD INSTALL --preclean --no-multiarch --no-test-load .")
```
to compile and install the functions.

## Usage

### ``simFM()`` -- Simulating a dynamic factor model

#### Parameters

- ``T`` -- Number of observations
- ``N`` -- Number of variables
- ``R`` -- Number of factors
- ``Sigma_epsilon`` -- Variance-covariance matrix of the transition errors
- ``Lambda`` -- Factor loadings matrix
- ``mu_xi`` -- Mean of the measurement error
- ``Sigma_xi`` -- Variance-covariance matrix of the measurement error
- ``Phi`` -- Factor VAR coefficient matrix (either a list of length $P$ of $(R\times R)$ matrices or a matrix of dimensions $(R\times RP)$)
- ``P`` -- Order of the factor VAR process
- ``quarterfy`` -- Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
- ``corr`` -- Indicating whether or not the measurement error should be internatlly, randomly, cross-crossectionally correlated (default is ``FALSE``)
- ``beta_param`` -- Beta parameter governing the degree of correlation of the measurement error (default is ``Inf``, i.e., disabled)
- ``m`` -- Ratio of monthly predictors ought to be quarterfied (default is ``0``)
- ``seed`` -- Seed for everything random (default is ``20022024``)
- ``burn_in`` -- Burn-in period (default is ``1000``)
- ``rescale`` -- Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component (default is ``TRUE``)
- ``check_stationarity`` -- Indicating whether the factor VAR process should be checked for stationarity (default is ``FALSE``)
- ``stationarity_check_threshold`` -- Threshold for checking for a unit root (default is ``1e-5``)

#### Example

```R

# Draw from an approximate mixed frequency factor model
T <- 600 # Number of observations
N <- 100 # Number of variabes
R <- 2 # Number of factors
Sigma_epsilon <- diag(1, R) # Variance-covariance matrix of the transition errors
Lambda <- matrix(rnorm(N * R), N, R) # Factor loadings matrix
mu_xi <- rep(0, N) # Mean of the measurement error
Sigma_xi <- diag(1, N) # Variance-covariance matrix of the measurement error
Phi <- cbind(diag(0.5, R), -diag(0.25, R)) # Factor VAR coefficient matrix
P <- 2 # Order of the factor VAR process
quarterfy <- TRUE # Indicating whether or not some of the variables should be aggregated to quarterly observations (i.e., quarterfied)
corr <- TRUE # Indicating whether or not the measurement error should be internatlly cross-crossectionally correlated
beta_param <- 1 # Beta parameter governing the degree of correlation of the measurement error
m <- 0.03 # Ratio of monthly predictors ought to be quarterfied
seed <- 16022024 # Seed
set.seed(seed)
burn_in <- 999 # Burn-in period
rescale <- TRUE # Indicating whether the variance of the measurement error should be scaled according to the variance of the common-component 

# Draw data
FM_mf <- simFM(T = T, N = N, R = R, Lambda = Lambda, mu_xi = mu_xi, Sigma_xi = Sigma_xi,
            Sigma_epsilon = Sigma_epsilon, Phi = Phi, P = P, quarterfy = quarterfy, m = m,
            corr = corr, beta_param = beta_param, seed = seed, burn_in = burn_in, rescale = rescale)


```

### ``TwoStepSDFM()`` -- Estimating the model parameters and higher-frequency factors

#### Parameters

- ``X`` -- $(N\times T)$ matrix of observations
- ``delay`` -- $(N \times 1)$ vector of publication lag (in months)
- ``selected`` -- $(R\times 1)$ vector of the number of non-zero loadings for each factor
- ``R`` -- Number of factors
- ``P`` -- Maximum lag length of the factor VAR processes (default is ``10``)
- ``decorr_errors`` -- Indicating whether the measurement error should be treated as cross-sectionally correlated (default is ``TRUE``)
- ``crit`` -- Information criterion to choose the factor VAR process lag length (available are ``"BIC"``, ``"AIC"``, ``"HC"``; default is ``"BIC"``)
- ``l2`` -- Ridge penalty (default is ``1e-6``)
- ``l1`` -- $(R\times 1)$ vector of LASSO penalties for each factor when using LARS (default is ``NaN``,i.e., disabled and ``selected`` is used as stopping criterion)
- ``max_iterations`` -- Maximum number of iterations (default is ``1000``)
- ``steps`` -- Maximum number of steps when using the LARS algorithm (default is ``NaN``, i.e., disabled and ``selected`` is used as stopping criterion)
- ``comp_null`` -- Computational zero (default is ``1e-15``)
- ``check_rank`` -- Indicating whether to check the rank of the variance-covariance matrix of the data (default is ``FALSE``)
- ``conv_crit`` -- Conversion criterion for the iterative SPCA algorithm (default is ``1e-4``)
- ``conv_threshold`` -- Conversion threshold when using coordinate descent (default is ``1e-4``)
- ``log`` -- Indicating whether to print output for monitoring purposes (default is ``FALSE``)

#### Example

```R

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
selected <- c(round(N * 0.8), round(N * 0.75))
fit_sparse <- twoStepSDFM(FM$X, delay, selected, R, l2 = 1e+4)

```

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

© 2024 Domenic Franjic

This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

This work is partially based on the LARS-EN and SPCA algorithms found in @Zou2020elasticnet.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

**To Contribute:**

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes with descriptive messages.
4. Push to your fork and submit a pull request.

## Support

If you have any questions or need assistance, please open an issue on the GitHub repository or contact us via email.

## Contact

- **Name**: Domenic Franjic
- **Institution**: University of Hohenheim
- **Department**: Econometrics and Statistics, Core Facility Hohenheim
- **E-Mail**: franjic@uni-hohenheim.de

### References

- Bates, Douglas, and Dirk Eddelbuettel. 2013. “Fast and Elegant Numerical Linear Algebra Using the RcppEigen Package.” Journal of Statistical Software 52 (5): 1–24.
- Doz, Catherine, Domenico Giannone, and Lucrezia Reichlin. 2011. “A two-step estimator for large approximate dynamic factor models based on Kalman filtering.” Journal of Econometrics 164 (1): 188–205.
- Eddelbuettel, Dirk, and Romain François. 2011. “Rcpp: Seamless R and C++ Integration.” Journal of Statistical Software 40 (8): 1–18.
- Franjic, Domenic, and Karsten Schweikert. 2024. “Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model.” Available at SSRN.
- Giannone, Domenico, Lucrezia Reichlin, and David Small. 2008. “Nowcasting: The Real-Time Informational Content of Macroeconomic Data.” Journal of Monetary Economics 55 (4): 665–76.
- Guennebaud, Gaël, Benoît Jacob, et al. 2010. “Eigen V3.” http://eigen.tuxfamily.org.
- Koopman, Siem Jan, and James Durbin. 2000. “Fast Filtering and Smoothing for Multivariate State Space Models.” Journal of Time Series Analysis 21 (3): 281–96.
- Mariano, Roberto S., and Yasutomo Murasawa. 2003. “A New Coincident Index of Business Cycles Based on Monthly and Quarterly Series.” Journal of Applied Econometrics 18 (4): 427–43.
- Zou, Hui, and Trevor Hastie. 2020. Elasticnet: Elastic-Net for Sparse Estimation and Sparse PCA. https://CRAN.R-project.org/package=elasticnet.
- Zou, Hui, Trevor Hastie, and Robert Tibshirani. 2006. “Sparse Principal Component Analysis.” Journal of Computational and Graphical Statistics 15 (2): 265–86.

