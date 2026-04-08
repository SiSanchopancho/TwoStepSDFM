#' Example factor model dataset
#'
#' This is a simulated factor model dataset for demonstration and testing
#' with `TwoStepSDFM`.
#'
#' @format A `SimulData` containing the following elements:
#' \describe{
#' \item{data}{If `starting_date` is provided, a `zoo` object, else, 
#' a (`no_of_vars` \eqn{\times}{x} `no_of_obs`) numeric matrix holding the 
#' simulated data.}
#' \item{factors}{If `starting_date` is provided, a `zoo` object, else a
#'  (`no_of_factors` \eqn{\times}{x} `no_of_obs`) numeric matrix holding  the 
#'  simulated latent factors.}
#' \item{trans_var_coeff}{Numeric (`no_of_factors` \eqn{\times}{x} 
#' (`no_of_factors` * `factor_lag_order`)) factor VAR coefficient matrix.}
#' \item{loading_matrix}{Numeric factor loading matrix.}
#' \item{meas_error}{If `starting_date` is provided, a `zoo` object, else 
#' a (`no_of_vars` \eqn{\times}{x} `no_of_obs`) numeric matrix holding the 
#' fundamental measurement errors.}
#' \item{meas_error_var_cov}{Numeric measurement error 
#' variance-covariance matrix.}
#' \item{trans_error_var_cov}{Numeric transition error 
#' variance-covariance matrix.}
#' \item{frequency}{Integer vector of variable frequencies.}
#' \item{delay}{Integer vector of variable delays, measured as the number of 
#' months since the latest available observation.}
#' }
#' 
#' @source Generated via `simFM()`. For details see `factor_model$call`.
"factor_model"

#' Mixed-frequency factor model dataset
#'
#' This dataset contains simulated mixed-frequency factor model data
#' for examples in `TwoStepSDFM`.
#'
#' @format A `SimulData` containing the following elements:
#' \describe{
#' \item{data}{If `starting_date` is provided, a `zoo` object, else, 
#' a (`no_of_vars` \eqn{\times}{x} `no_of_obs`) numeric matrix holding the 
#' simulated data.}
#' \item{factors}{If `starting_date` is provided, a `zoo` object, else a
#'  (`no_of_factors` \eqn{\times}{x} `no_of_obs`) numeric matrix holding  the 
#'  simulated latent factors.}
#' \item{trans_var_coeff}{Numeric (`no_of_factors` \eqn{\times}{x} 
#' (`no_of_factors` * `factor_lag_order`)) factor VAR coefficient matrix.}
#' \item{loading_matrix}{Numeric factor loading matrix.}
#' \item{meas_error}{If `starting_date` is provided, a `zoo` object, else 
#' a (`no_of_vars` \eqn{\times}{x} `no_of_obs`) numeric matrix holding the 
#' fundamental measurement errors.}
#' \item{meas_error_var_cov}{Numeric measurement error 
#' variance-covariance matrix.}
#' \item{trans_error_var_cov}{Numeric transition error 
#' variance-covariance matrix.}
#' \item{frequency}{Integer vector of variable frequencies.}
#' \item{delay}{Integer vector of variable delays, measured as the number of 
#' months since the latest available observation.}
#' }
#' 
#' @source Generated via Generated via `simFM()`. For details see 
#' `mixed_freq_factor_model$call`.
"mixed_freq_factor_model"