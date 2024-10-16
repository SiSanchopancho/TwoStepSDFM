# TwoStepSDFM
A ``C++``-based ``R`` implementation of a two-step estimation procedure for a (linear Gaussian) Sparse Dynamic Factor Model (SDFM) as outlined in Franjic and Schweikert (2024).

## Introduction

The ``TwoStepSDFM`` package provides a fast implementation of the Kalman Filter and Smoother (KFS) to estimate factors in a mixed-frequency SDFM framework, explicitly accounting for cross-sectional correlation in the measurement error. The KFS is initialized using results from Sparse Principal Components Analysis (SPCA) by Zou and Hastie (2006) in a preliminary step. This approach generalizes the two-step methods by Giannone et al. (2008) and Doz et al. (2011) by allowing for a potentially dense factor loadings matrix.

We estimate the loadings matrix Λ and factors fₜ of the following state-space model:

xₜ = Λ fₜ + ξₜ
fₜ = Σₚ Φₚ fₜ₋ₚ + εₜ

for t = 1, ..., T. We use SPCA to estimate Λ, resulting in exact zero elements due to regularization. We also explicitly model the non-diagonal variance-covariance matrix of the measurement error ξₜ.

