# TwoStepSDFM
A ``C++``-based ``R`` implementation of a two-step estimation procedure for a (linear Gaussian) Sparse Dynamic Factor Model (SDFM) as outlined in Franjic and Schweikert (2024).

## Introduction

The ``TwoStepSDFM`` package provides a fast implementation of the Kalman Filter and Smoother (KFS) to estimate factors in a mixed-frequency SDFM framework, explicitly accounting for cross-sectional correlation in the measurement error. The KFS is initialized using results from Sparse Principal Components Analysis (SPCA) by Zou and Hastie (2006) in a preliminary step. This approach generalizes the two-step estimator for approximate dynamic factor models by Giannone, Reichlin, and Small (2008) and Doz, Giannone, and Reichlin (2011).
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

