#' ## usethis namespace: start
#' @importFrom Rcpp sourceCpp
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


#' Run Sparse DFM KFS
#'
#' This function wraps the C++ implementation of the Sparse DFM KFS to
#' provide default parameters and potentially pre-process inputs or
#' post-process outputs.
#'
#' @param X Numeric "N x T" matrix of data.
#' @param delay Integer vector of predictor delays.
#' @param selected Integer vector of the number of selected variables for each factor.
#' @param R Integer, number of factors.
#' @param P Integer, max P of the VAR(P) process of the factors.
#' @param decorr_errors Logical, whether or not the errors should be decorrelated.
#' @param crit Character string, information criterion used for the estimation.
#' @param method Character string, method used for estimation ("BIC", "AIC", "HIC").
#' @param l2 Numeric, ridge penalty.
#' @param l1 Numeric vector, lasso penalties for each factor (set to NaN if not intended as stopping criterion).
#' @param max_iterations Integer, maximum number of iterations.
#' @param steps Integer, number of steps (used for LARS-EN as an alternative).
#' @param comp_null Numeric, computational zero.
#' @param check_rank Logical, whether or not the rank of the variance-covariance matrix should be checked.
#' @param conv_crit Numeric, conversion criterion for the SPCA algorithm.
#' @param conv_threshold Numeric, conversion criterion for the coordinate descent algorithm.
#' @param log Logical, whether or not output should be printed along the algorithm.
#' @return 
#' The `twoStepSDFM` function returns an object `fit` which, contains the following elements:
#' \describe{
#'   \item{\code{Lambda_hat}}{A matrix of estimated loadings for each factor on each observed variable.}
#'   \item{\code{F}}{The factor estimates.}
#'   \item{\code{Pt}}{The estimated precision matrix of the factors.}
#'   \item{\code{Wt}}{The estimated variance-covariance matrix of the residuals/errors.}
#'   \item{\code{C}}{A matrix representing any transformations applied to the factors for error decorrelation, if \code{decorr_errors} is TRUE.}
#'   \item{\code{P}}{The order of the VAR(P) model used in the factor model process.}
#' }
#' @export
#' @examples
#' \dontrun{
#' result <- TwoStepSDFM(X = matrix_data, delay = delays, selected = selections)
#' }
twoStepSDFM <- function(X,
                        delay,
                        selected,
                        R,
                        P = 10,
                        decorr_errors = TRUE,
                        crit = "BIC",
                        l2 = 1e-6,
                        l1 = NaN,
                        max_iterations = 1000,
                        steps = NaN,
                        comp_null = 1e-15,
                        check_rank = FALSE,
                        conv_crit = 1e-4,
                        conv_threshold = 1e-4,
                        log = FALSE) {
  # Misshandling of the loadings matrix
  
  XR <- try(as.matrix(X), silent = TRUE)
  if (inherits(XR, "try-error")) {
    stop(paste0("X must be a matrix or convertable to a matrix"))
  }
  
  XR <- t(XR)
  
  # Misshandling of selected
  
  if (length(selected) == 1) {
    if (is.nan(selected) || is.null(selected) || is.infinite(-selected)) {
      selected <- rep(dim(XR)[2], R)
    }
  }
  
  if (length(selected) != R || any(selected < 0)) {
    stop(paste0("The vector containing the number of variables ought to be selected is either not of length R = ", R, " or has negative values. If you aim to disable selected as stopping criterion set it to NaN, -Inf, or Null."))
  }
  
  l1_disabled <- is.null(l1) || (length(l1) == 1 && (is.nan(l1) || is.infinite(-l1)))
  if (!l1_disabled && (length(l1) != R || any(l1 < 0))) {
    stop(paste0("The vector of l1 penalties is either not of length R = ", R, " or has negative values. If you aim to disable the l1 penalty as stopping criterin set it to NaN, -Inf, or Null."))
  }
  
  if (l1_disabled) {
    l1 <- rep(-2147483647L, R)
  }
  
  if (!is.na(steps) && !is.null(steps)) {
    if (steps < 0 && !is.infinite(-steps)) {
      stop(paste0("The maximum number of steps used for the LARS-EN algorith must be non-negative. To disable the number of steps as stopping criterion, set \"steps\" to NaN, -Inf, or Null"))
    }
  }
  
  if (is.nan(steps) || is.null(steps) || is.infinite(-steps)) {
    steps <- -2147483647L # C++ INT_MIN
  }
  
  if (steps != floor(steps) && !is.na(steps)) {
    stop(paste0("The maximum number of steps used for the LARS-EN algorith must be an integer."))
  }
  
  # Disable the checking of the KFS converions
  
  KFS_conv_crit <- 1e+15
  
  # Call the C++ function
  
  result <- runSDFMKFS(
    X_in = XR,
    delay = delay,
    selected = selected,
    R = R,
    order = P,
    decorr_errors = decorr_errors,
    crit = crit,
    l2 = l2,
    l1 = l1,
    max_iterations = max_iterations,
    steps = steps,
    comp_null = comp_null,
    check_rank = check_rank,
    conv_crit = conv_crit,
    conv_threshold = conv_threshold,
    log = log,
    KFS_conv_crit = KFS_conv_crit
  )
  
  # Return the result
  return(result)
}


#' Draw data from a dynamic factor model
#'
#' The function simulates an (approximate) Dynamic Factor Model
#'
#' @param T Number of observations.
#' @param N Number of Variables.
#' @param R Number of factors.
#' @param Lambda Loadings matrix.
#' @param mu_xi Mean of the measurement errors.
#' @param Sigma_xi Variance-covariance matrix of the measurement errors.
#' @param Sigma_epsilon Variance-covariance matrix of the transition errors.
#' @param Phi Either a list of length P or a matrix of dimensions Rx(R * P) or (R * P)xR holding the VAR coefficients of the factor VAR process.
#' @param P Order of the VAR process in the state equation.
#' @param quarterfy Whether or not some of the data should be aggregated to quarterly representations.
#' @param m Ratio of variables ought to be quarterfied.
#' @param corr Whether or not the measurement error should be randomly correlated inside the function using a random correlation matrix with off-diagonal elements governed by a beta-distribution.
#' @param beta_param Parameter of the beta-distribution governing the off-diagonal elements of the variance-covariance matrix of the measurement error.
#' @param seed Seed for all random processes inside the function.
#' @param burn_in Burn-in period of the simulated data ought to be discarded at the beginning of the data.
#' @param rescale Whether or not the variance of the measurement error should be rescaled by the common component to equalise the signal-to-noise ratio.
#' @param check_stationarity Whether or not the stationarity properties of the factor VAR process should be checked.
#' @param stationarity_check_threshold Threshold of the stationarity check for when to deem an eigenvalue negative beyond the numerical error.
#' @return Returns an `simFM` object containing the DFM parameters expressed in the factor VAR companion form.
#' The object includes:
#' \describe{
#'   \item{\code{F}}{A \eqn{RP \times T} matrix of simulated factors.}
#'   \item{\code{Phi}}{The \eqn{R \times RP} factor VAR coefficients matrix.}
#'   \item{\code{Lambda}}{The \eqn{N \times RP} loadings matrix.}
#'   \item{\code{Sigma_xi}}{The \eqn{N \times N} variance-covariance matrix of the measurement errors.}
#'   \item{\code{Sigma_epsilon}}{The \eqn{R \times R} variance-covariance matrix of the transition errors.}
#'   \item{\code{Xi}}{A \eqn{N \times T} matrix representing the measurement errors.}
#'   \item{\code{X}}{The \eqn{N \times T} matrix of simulated observed variables.}
#'   \item{\code{frequency}}{The \eqn{N \times 1} vector of frequency of the data.}
#' }
#' @export
simFM <- function(T, N, R, Lambda, mu_xi, Sigma_xi, Sigma_epsilon, Phi, P, quarterfy, m = 0, corr = FALSE,
                  beta_param = Inf, seed = 20022024, burn_in = 1000, rescale = TRUE,
                  check_staionarity = FALSE, stationarity_check_threshold = 1e-5) {
  # Check for miss-handling
  
  if(!is.logical(quarterfy) && !(quarterfy %in% c(0, 1))){
    
    stop(paste0("quarterfy should be logical."))
    
  }
  
  if(!is.logical(rescale) && !(rescale %in% c(0, 1))){
    
    stop(paste0("rescale should be logical."))
    
  }
  
  # Misshandling of the State-Space Dimensions
  
  if (!is.numeric(N) || !is.numeric(T) || !is.numeric(R) || !is.numeric(P)) {
    stop(paste0("Either the number of variables N, the number of Observations T, the number of Factors R, or the facgtor VAR order P is not provided as numeric."))
  }
  
  if (N <= 0 || T <= 0 || R <= 0 || P <= 0) {
    stop(paste0("Either the number of variables N, the number of Observations T, the number of Factors R, or the factor VAR process order P is negative or zero."))
  }
  
  if (max(c(N, T)) <= R) {
    warning(paste0("R = ", R, " >= ", max(N, T), " = max(N,T)! In general, R << N,T."))
  }
  
  # Misshandling of the loadings matrix
  
  LambdaR <- try(as.matrix(Lambda), silent = TRUE)
  if (inherits(LambdaR, "try-error")) {
    stop(paste0("Lambda must be a matrix or convertable to a matrix"))
  }
  
  if (dim(LambdaR)[1] != N || dim(LambdaR)[2] != R) {
    stop(paste0("The loadings Matrix Lambda must be of dimensions NxR = ", N, "x", R, " but is ", dim(LambdaR)[1], "x", dim(LambdaR)[2], "."))
  }
  
  # Misshandling of the measurement error mean
  
  muR <- try(as.matrix(mu_xi), silent = TRUE)
  if (inherits(muR, "try-error")) {
    stop(paste0("mu_xi must be a vector or convertable to a matrix"))
  }
  
  if (dim(muR)[1] == 1 && dim(muR)[2] == N) {
    muR <- t(muR)
  }
  
  if(!is.numeric(muR)){
    
    stop(paste0("mu_xi must be of type numeric"))
    
  }
  
  if (dim(muR)[1] != N && dim(muR)[2] != 1) {
    stop(paste0("The mean of the measurement error mu_xi must be a vector of dimensions ", N, "x1 or 1x", N, " but is ", dim(muR)[1], "x", dim(muR)[2], "."))
  }
  
  # Misshandling of the measurement error variance-covariance matrix
  
  SigmaR <- try(as.matrix(Sigma_xi), silent = TRUE)
  if (inherits(SigmaR, "try-error")) {
    stop(paste0("Sigma_xi must be a matrix or convertable to a matrix"))
  }
  
  if(!is.numeric(SigmaR)){
    
    stop(paste0("Sigma_xi must be of type numeric"))
    
  }
  
  if (sum(dim(SigmaR) == c(N, N)) != 2) {
    stop(paste0("The variance-covariance matrix of the measurement error Sigma_xi must be a vector of dimensions ", N, "x", N, " but is ", dim(SigmaR)[1], "x", dim(SigmaR)[2], "."))
  }
  
  if (!isSymmetric(SigmaR)) {
    stop(paste0("The variance-covariance matrix of the measurement error must be symmetric."))
  }
  
  # Misshandling of the transition error variance-covariance matrix
  
  SR <- try(as.matrix(Sigma_epsilon), silent = TRUE)
  if (inherits(SR, "try-error")) {
    stop(paste0("Sigma_epsilon must be a matrix or convertable to a matrix"))
  }
  
  if (!all(dim(SR) == c(R, R))) {
    stop(paste0("The variance-covariance matrix of the transition error Sigma_epsilon must be a vector of dimensions ", R, "x", R, " but is ", dim(SigmaR)[1], "x", dim(SigmaR)[2], "."))
  }
  
  if (!isSymmetric(SR)) {
    stop(paste0("The variance-covariance matrix of the transition error must be symmetric."))
  }
  
  if (!isSymmetric(SR)) {
    stop(paste0("The variance-covariance matrix of the transition error must be symmetric."))
  }
  
  # Misshandling of the VAR coefficient matrix
  
  if (is.list(Phi) && !is.data.frame(Phi)) {
    # If Phi is provided as list
    
    # Check whether each element in Phi has the correct dimensions and whether there is the correct number of matrices provided
    
    ind <- c()
    for (i in 1:length(Phi)) {
      s <- sum(dim(Phi[[i]]) == c(R, R))
      if (s != 2) {
        ind[i] <- i
      }
    }
    
    if (!is.null(ind)) {
      stop("The VAR coefficient matrices in Phi must be of dimensions ", R, "x", R, ", which is not the case for the matrices ", paste(ind, collapse = ", "), ".")
    }
    
    if (length(Phi) != P) {
      stop("The number of VAR coefficient matrices in Phi must equal P = ", P, " but is ", length(Phi), ".")
    }
    
    # Store the list elements as matrix
    
    PhiR <- matrix(0, R, R * P)
    for (o in 1:P) {
      PhiR[1:R, ((o - 1) * R + 1):((o - 1) * R + R)] <- Phi[[o]]
    }
  } else {
    # If Phi is provided as matrix or equivalent
    
    PhiR <- try(as.matrix(Phi), silent = TRUE)
    if (inherits(PhiR, "try-error")) {
      stop(paste0("PhiR must be a list, a matrix or convertable to a matrix"))
    }
    
    if (all(dim(PhiR) == c(P * R, R))) {
      PhiR <- t(PhiR)
    }
    
    if (!all(dim(PhiR) == c(R, R * P))) {
      stop(paste0("The VAR coefficient matrix of the factor process must be a vector of dimensions ", R, "x", R * P, " or ", R * P, "x", R, ", where each ", R, "x", R, " block represents a VAR coefficient of corresponding lag P. Here, the provided matrix is ", dim(PhiR)[1], "x", dim(PhiR)[2], "."))
    }
  }
  
  if(!is.numeric(PhiR)){
    
    stop("Phi must be a numeric matrix or a list of numeric matrices")
    
  }
  
  # Check whether the process will result in a stationary process
  
  if (check_staionarity) {
    
    Comp <- matrix(0, R * P, R * P)
    
    Comp[1:R, 1:(P * R)] <- PhiR
    
    Comp[(R + 1):(P * R), 1:(R * (P - 1))] <- diag(1, R * (P - 1))
    
    # Compute eigenvalues
    e <- eigen(Comp)$values
    test <- any(abs(abs(e) - 1) < stationarity_check_threshold)
    
    if (test) {
      warning("At least one eigenvalue very close to one detected. Process might have random walk properties.\n")
    }
    
    test2 <- all(abs(e) < 1 + stationarity_check_threshold)
    
    if (!test2) {
      warning("At least one eigenvalue lies outside the complex unit-circle. Process is not stationary.")
    }
    
  }
  
  # Misshandling of the beta parameter
  
  if (beta_param <= 0) {
    stop("The value of the parameter of the beta distribution governing the cross-sectional correlation in the measurement error beta_param must be strictly positive. To disable the random cross-section correlation of the measurement error set corr = FALSE.")
  }
  
  if (corr && is.infinite(beta_param)) {
    stop("If the data should be randomly correlated internally beta_param should not be Inf. To disable the random cross-section correlation of the measurement error set corr = FALSE.")
  }
  
  # M;isshandling of the ratio of quarterly variables
  
  if (!is.numeric(m)) {
    stop("The ratio of variables ought to be quarterfied m must be numeric.")
  }
  
  if (m < 0) {
    stop("The ratio of variables ought to be quarterfied m must be non-negative.")
  }
  
  # Misshandling of the burn_in period
  
  if (burn_in < 0) {
    stop("The burn in period at the start of the simulated data burn_in must be non-negative.")
  }
  
  FM <- runStaticFM(
    T = T, N = N, S = SR, Lambda = LambdaR, mu_e = muR, Sigma_e = SigmaR, A = PhiR, order = P, quarterfy = quarterfy, corr = corr, beta_param = beta_param,
    m = m, seed = seed, R = R, burn_in = burn_in, rescale = rescale
  )
  
  return(FM)
}
