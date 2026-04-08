/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright (C) 2024 Domenic Franjic
 *
 * This file is part of TwoStepSDFM.
 *
 * TwoStepSDFM is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * TwoStepSDFM is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.
 */

#define _USE_MATH_DEFINES // If you need some math constants

 // Externakl includes
#include <random>
#include <Eigen/Eigen>
#include <math.h>
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "TwoStepSDFM_types.h"
#include "Internals/DataGen.h"
#include "Internals/SparseDFM.h"


//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `twoStepSDFM()` instead for a stable and
//' supported interface.
//'
// [[Rcpp::export]]
Rcpp::List runUVKFS(
  Rcpp::NumericMatrix X_in,
  Rcpp::IntegerVector delay,
  Rcpp::NumericMatrix state_var_cov,
  Rcpp::NumericMatrix measurement_var_cov,
  Rcpp::NumericMatrix loading_matrix,
  Rcpp::NumericMatrix factor_var_coefficient_matrices,
  const int R,
  const int order,
  const int fcast_horizon,
  const bool decorr_errors,
  double comp_null,
  const bool parallel,
  const double jitter
)
{

  // Map the numeric matrices and vectors to eigen objects
  Eigen::Map<Eigen::MatrixXd> X_in_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X_in));
  Eigen::Map<Eigen::VectorXi> delay_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delay));
  Eigen::Map<Eigen::MatrixXd> state_var_cov_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(state_var_cov));
  Eigen::Map<Eigen::MatrixXd> measurement_var_cov_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(measurement_var_cov));
  Eigen::Map<Eigen::MatrixXd> loading_matrix_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(loading_matrix));
  Eigen::Map<Eigen::MatrixXd> factor_var_coefficient_matrices_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(factor_var_coefficient_matrices));

  // Enable/disable parallelisation in Eigen
  if (parallel) {
    Eigen::setNbThreads(0);
  }
  else {
    Eigen::setNbThreads(1);
  }

  SparseDFM::SDFM<SparseDFM::Structure::SPARSE> results(X_in_eigen, R, order);
  results.loading_matrix = loading_matrix_eigen;

  // Build companion form matrices and cut data to size
  Eigen::MatrixXd effective_data = results.data.transpose();
  int ind_of_first_companion_obs = results.max_factor_var_order - 1;
  int no_of_states = results.no_of_factors * results.max_factor_var_order;
  Eigen::MatrixXd filter_var_cov_initial = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
  filter_var_cov_initial.diagonal().setConstant(1000);
  Eigen::MatrixXd comp_form_loading_matrix = Eigen::MatrixXd::Zero(results.no_of_vars, no_of_states);
  Eigen::MatrixXd state_equation_var_cov = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
  Eigen::MatrixXd meas_equation_var_cov = Eigen::MatrixXd::Zero(results.no_of_vars, results.no_of_vars);
  Eigen::MatrixXd comp_form_factor_var_coeff = Eigen::MatrixXd::Zero(no_of_states, no_of_states);

  // Estimate the missing state and measurement equation parameters
  int index_of_full_sample = results.no_of_obs - delay_eigen.maxCoeff() - 1;
  {

    // State equation parameters
    Eigen::MatrixXd ident_no_of_states = Eigen::MatrixXd::Identity(no_of_states, no_of_states);
    comp_form_factor_var_coeff.topLeftCorner(results.no_of_factors, results.no_of_factors * results.max_factor_var_order) = factor_var_coefficient_matrices_eigen;
    comp_form_factor_var_coeff.bottomLeftCorner(results.no_of_factors * (ind_of_first_companion_obs), results.no_of_factors * (ind_of_first_companion_obs)).setIdentity();
    comp_form_factor_var_coeff.bottomRightCorner(results.no_of_factors * (ind_of_first_companion_obs), results.no_of_factors).setZero();
    Eigen::MatrixXd factor_var_coeff_transposed = comp_form_factor_var_coeff.transpose();
    state_equation_var_cov.topLeftCorner(results.no_of_factors, results.no_of_factors) = state_var_cov_eigen;

    // Measurement equation parameters
    comp_form_loading_matrix.topLeftCorner(results.no_of_vars, results.no_of_factors) = results.loading_matrix;
    meas_equation_var_cov = measurement_var_cov_eigen;
    meas_equation_var_cov.diagonal().array() += jitter; // jitter the matrix for numerical stability
  }
  results.factors.conservativeResize(no_of_states, (results.no_of_obs + fcast_horizon + 1));
  results.factors.setZero();
  results.filter_var_cov.conservativeResize(no_of_states, no_of_states * (results.no_of_obs + fcast_horizon + 1));
  results.filter_var_cov.setZero();
  results.smoother_var_cov.conservativeResize(no_of_states, no_of_states * results.no_of_obs);
  results.smoother_var_cov.setZero();
  results.zero_indeces.conservativeResize(results.loading_matrix.rows(), results.loading_matrix.cols());
  results.zero_indeces.setZero();

  if (decorr_errors) { // Decorrelate the data and filter dynamics for the implementation of the univariate Kalman filter and smoother
    results.zero_indeces = results.loading_matrix.unaryExpr([comp_null](double a) { return (a == 0) ? 1 : 0; });
    Eigen::LLT<Eigen::MatrixXd> meas_llt(meas_equation_var_cov);
    if (meas_llt.info() != Eigen::Success) {
#if !defined(_MSC_VER)
      Rcpp::Rcerr << "\nWARNING: Decorrelation failed. Using robust cholesky decomposition.\n";
#else
      std::cerr << "\nWARNING: Decorrelation failed. Using robust cholesky decomposition.\n";
#endif
      Eigen::LDLT<Eigen::MatrixXd> meas_ldlt(meas_equation_var_cov);
      if (meas_ldlt.info() != Eigen::Success) {
#if !defined(_MSC_VER)
        Rcpp::Rcerr << "\nWARNING: Robust decorrelation failed. Ignoring measurement error correlation structure.\n";
#else
        std::cerr << "\nWARNING: Robust decorrelation failed. Ignoring measurement error correlation structure.\n";
#endif
        results.inv_chol_variable_var_cov = Eigen::MatrixXd::Identity(results.no_of_vars, results.no_of_vars);
      }
      Eigen::MatrixXd lower_factor = meas_ldlt.matrixL();
      Eigen::VectorXd diagonal_matrix_diag = meas_ldlt.vectorD();
      diagonal_matrix_diag.diagonal() = diagonal_matrix_diag.diagonal().unaryExpr([comp_null](double x) { return std::max(x, comp_null); }).eval();
      Eigen::DiagonalMatrix<double, Eigen::Dynamic> sqrt_diagonal_matrix_diag = meas_ldlt.vectorD().cwiseSqrt().asDiagonal();
      Eigen::Transpositions<Eigen::Dynamic> permutation_matrix = meas_ldlt.transpositionsP();
      Eigen::MatrixXd chol_variable_var_cov = lower_factor * sqrt_diagonal_matrix_diag;
      results.inv_chol_variable_var_cov = chol_variable_var_cov.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(results.no_of_vars, results.no_of_vars));;
    }
    else {
      results.inv_chol_variable_var_cov = meas_llt.matrixL().solve(Eigen::MatrixXd::Identity(results.no_of_vars, results.no_of_vars));
    }

    if (0 < delay_eigen.maxCoeff()) { // Decorrelate row by row in the case of missing data
      Eigen::MatrixXd Temp = effective_data;
      for (int curr_obs = 0; curr_obs < results.no_of_obs; ++curr_obs) {
        Eigen::ArrayXd missing_ind = (!((Eigen::VectorXi::Constant(results.no_of_vars, 1, results.no_of_obs) - delay_eigen).array() <= curr_obs)).cast<double>();
        for (int curr_var = 0; curr_var < results.no_of_vars; ++curr_var)
        {
          if (results.no_of_obs - delay_eigen(curr_var) <= curr_obs) {
            Temp(curr_var, curr_obs) = DBL_MAX;
          }
          else {
            Eigen::VectorXd Cn = (results.inv_chol_variable_var_cov.row(curr_var).transpose().array() * missing_ind).matrix();
            Eigen::VectorXd Xt = effective_data.col(curr_obs).array().isNaN().select(0.0, effective_data.col(curr_obs)).matrix();
            Temp(curr_var, curr_obs) = Cn.dot(Xt);
          }
        }
      }
      effective_data = Temp;
    }
    else
    {
      effective_data = results.inv_chol_variable_var_cov * effective_data;
    }
    comp_form_loading_matrix.topLeftCorner(results.no_of_vars, results.no_of_factors) = (results.inv_chol_variable_var_cov * comp_form_loading_matrix.topLeftCorner(results.no_of_vars, results.no_of_factors)).eval();
    meas_equation_var_cov = (results.inv_chol_variable_var_cov * meas_equation_var_cov * results.inv_chol_variable_var_cov.transpose()).eval();

  }
  else { // Store identity as inverse cholesky factor in the case of no decorrelation
    results.inv_chol_variable_var_cov = Eigen::MatrixXd::Identity(results.no_of_vars, results.no_of_vars);
  }

  /* State filtering and smoothing */

  Eigen::MatrixXd inv_meas_equation_var_cov = Eigen::MatrixXd::Zero(results.no_of_vars, results.no_of_obs);
  Eigen::MatrixXd kalman_gain = Eigen::MatrixXd::Zero(no_of_states, results.no_of_vars * results.no_of_obs);
  Eigen::MatrixXd filter_error = Eigen::MatrixXd::Zero(results.no_of_vars, results.no_of_obs);

  int kalman_filter_starting_index = ind_of_first_companion_obs;
  Eigen::MatrixXd comp_form_factors = Eigen::MatrixXd::Zero(results.no_of_factors * results.max_factor_var_order, results.no_of_obs - delay_eigen.maxCoeff());
  results.univariateRepOfMultivariateKalmanFilter(comp_form_factors, inv_meas_equation_var_cov, kalman_gain, filter_error, kalman_filter_starting_index, results.no_of_obs, no_of_states, effective_data, comp_form_loading_matrix, comp_form_factor_var_coeff, meas_equation_var_cov, state_equation_var_cov, filter_var_cov_initial, delay_eigen, fcast_horizon);
  results.univariateRepOfMultivariateKalmanSmoother(kalman_gain, inv_meas_equation_var_cov, filter_error, results.no_of_obs, no_of_states, effective_data, comp_form_loading_matrix, comp_form_factor_var_coeff, meas_equation_var_cov, state_equation_var_cov, delay_eigen);
  results.loading_matrix = comp_form_loading_matrix;

  Eigen::setNbThreads(0);

  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("Pt") = Rcpp::wrap(results.filter_var_cov),
    Rcpp::Named("F") = Rcpp::wrap(results.factors),
    Rcpp::Named("Wt") = Rcpp::wrap(results.smoother_var_cov),
    Rcpp::Named("llt_success_code") = Rcpp::wrap(results.llt_success_code)
  );
}