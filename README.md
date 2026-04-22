# TwoStepSDFM
A ``C++``-based ``R`` implementation of the two-step estimation procedure for a (linear Gaussian) Sparse Dynamic Factor Model (SDFM) as outlined in Franjic and Schweikert (2024).

## Introduction

The ``TwoStepSDFM`` package provides a fast implementation of the Kalman Filter and Smoother (hereinafter KFS, see Koopman and Durbin, 2000) to estimate factors in a mixed-frequency SDFM framework, explicitly accounting for cross-sectional correlation in the measurement error. The KFS is initialized using results from Sparse Principal Components Analysis (SPCA) by Zou and Hastie (2006) in a preliminary step. This approach generalizes the two-step estimator for approximate dynamic factor models by Giannone, Reichlin, and Small (2008) and Doz, Giannone, and Reichlin (2011). For more details see Franjic and Schweikert (2024).

## Main Features

- **Fast Model Simulation**: The ``simFM()`` function provides a flexible framework to simulate mixed-frequency data with ragged edges from an approximate DFM.
- **Estimation of the Number of Factors**: The ``noOfFactors()`` function uses the Onatski (2009) procedure to estimate the number of factors efficiently while providing good finite sample performance.
- **Fast Model Estimation**: The ``twoStepSDFM()`` function provides a fast, memory-efficient, and convenient implementation of the two-step estimator outlined in Franjic and Schweikert (2024).
- **Fast Hyper-Parameter Cross-Validation**: The ``crossVal()`` function provides a fast and parallel cross-validation wrapper to retrieve the optimal hyper-parameters using time-series cross-validation (Hyndman and Athanasopoulos 2018) with random hyper-parameter search (Bergstra and Bengio 2012).
- **Fast Model Prediction**: The ``nowcast()`` function is a highly convenient prediction function for backcasts, nowcasts, and forecasts of multiple targets. It automatically takes care of all issues arising with mixed-frequency data and ragged edges.
- **Compatibility**: All functions take advantage of C++ for enhanced speed and memory-efficiency.

## Side Features

- **Fast dense DFM estimation and prediction**: The ``nowcast()`` function is also able to produce predictions of a dense DFM according to Giannone, Reichlin, and Small (2008). The function ``twoStepDenseDFM()`` additionally exposes an estimation procedure for the dense two-step estimator.
- **Fast SPCA**: ``sparsePCA()`` exposes the internal ``C++``-backed SPCA routine in ``R``. This provides access to a fast and memory-efficient SPCA estimation routine as implemented by Zou and Hastie (2020) in pure R.
- **Fast Kalman Filter and Smoother**: The ``kalmanFilterSmoother()`` function exposes the internal ``C++``-backed KFS routine.

## Installation

The package is available on CRAN and can be installed via ``install.packages("TwoStepSDFM")``. If this turns out to be no longer possible, run the PackageBuilder.R file.

## Prerequisites

For the installation from source via the PackageBuilder.R file, the following is required:

- **Rcpp**: A package for integrating `C++` code into `R` (Eddelbuettel and François, 2011). [Rcpp CRAN repository](https://CRAN.R-project.org/package=Rcpp)
- **RcppEigen**: A package for integrating the `Eigen` linear algebra library into `R` (Bates and Eddelbuettel, 2013). [RcppEigen CRAN repository](https://CRAN.R-project.org/package=RcppEigen)
- **GCC compiler** (version 5.0 or later) [GCC Website](https://gcc.gnu.org/).

## Usage

For a quick step-by-step user guide of the main features, see the package vignette.

## License

License: GPL v3

(C) 2024-2026 Domenic Franjic

This project is licensed under the **GNU General Public License v3.0**. See the LICENSE file for details.

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

## References

### Papers
- Bergstra, James, and Yoshua Bengio (2012). [Random Search for Hyper-Parameter Optimization](http://www.jmlr.org/papers/volume13/bergstra12a/bergstra12a.pdf). *Journal of Machine Learning Research*, 13(2).
- Doz, Catherine, Domenico Giannone, and Lucrezia Reichlin (2011). [A two-step estimator for large approximate dynamic factor models based on Kalman filtering](https://doi.org/10.1016/j.jeconom.2011.02.012). *Journal of Econometrics*, 164(1), 188–205.
- Franjic, Domenic, and Karsten Schweikert (2024). [Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model](https://ssrn.com/abstract=4733872). SSRN 4733872.
- Giannone, Domenico, Lucrezia Reichlin, and David Small (2008). [Nowcasting: The Real-Time Informational Content of Macroeconomic Data](https://doi.org/10.1016/j.jmoneco.2008.05.010). *Journal of Monetary Economics*, 55(4), 665–76.
- Koopman, Siem Jan, and James Durbin (2000). [Fast Filtering and Smoothing for Multivariate State Space Models](https://doi.org/10.1111/1467-9892.00186). *Journal of Time Series Analysis*, 21(3), 281–96.
- Marcellino, Massimiliano, and Christian Schumacher (2010). [Factor MIDAS for nowcasting and forecasting with ragged-edge data: A model comparison for German GDP](https://doi.org/10.1111/j.1468-0084.2010.00591.x). *Oxford Bulletin of Economics and Statistics*, 72(4), 518–50.
- Mariano, Roberto S., and Yasutomo Murasawa (2003). [A New Coincident Index of Business Cycles Based on Monthly and Quarterly Series](https://doi.org/10.1002/jae.695). *Journal of Applied Econometrics*, 18(4), 427–43.
- Onatski, Alexei (2009). [Testing Hypotheses about the Number of Factors in Large Factor Models](https://doi.org/10.3982/ECTA6964). *Econometrica*, 77(5), 1447–79.
- Zou, Hui, Trevor Hastie, and Robert Tibshirani (2006). [Sparse Principal Component Analysis](https://doi.org/10.1198/106186006X113430). *Journal of Computational and Graphical Statistics*, 15(2), 265–86.

### Books
- Hyndman, Rob J., and George Athanasopoulos (2018). *Forecasting: Principles and Practice* (3rd ed.). [OTexts Melbourne](https://otexts.com/fpp3/).

### Software / Packages
- Zou, Hui, and Trevor Hastie (2020). [Elasticnet: Elastic-Net for Sparse Estimation and Sparse PCA (R package)](https://CRAN.R-project.org/package=elasticnet)






