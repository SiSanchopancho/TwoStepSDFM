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

//void findLineOfProblem(int x, std::string additional_message = "") {
//  std::ofstream log_file;
//  log_file.open("debug_cpp.log", std::ios::app);
//  log_file << "Line executed: " << x << " " << additional_message << std::endl;
//  log_file.close();
//}
//
//template <typename Derived>
//void sendMatrixStraightToTheBank(const Eigen::MatrixBase<Derived>& matrix, std::string name = "matrix") {
//  std::ofstream log_file;
//  std::string file_name = name + ".log";
//  log_file.open(file_name, std::ios::app);
//  log_file << name << ": \n" << matrix << std::endl;
//  log_file.close();
//}

/* Function to run the two-step SDFM estimation procedure*/
/*Source:
 -  Franjic, Domenic and Schweikert, Karsten, Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model (February 21, 2024). Available at SSRN: https://ssrn.com/abstract=4733872 or http://dx.doi.org/10.2139/ssrn.4733872
 */

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `twoStepSDFM()` instead for a stable and
//' supported interface.
//'
// [[Rcpp::export]]
Rcpp::List runSDFMKFS(
  Rcpp::NumericMatrix X_in,
  Rcpp::IntegerVector delay,
  Rcpp::IntegerVector selected,
  int R,
  int order,
  bool decorr_errors,
  const char* crit,
  double l2,
  Rcpp::NumericVector l1,
  int max_iterations,
  int steps,
  Rcpp::NumericVector weights,
  double comp_null,
  double spca_conv_crit,
  const bool parallel,
  const unsigned fcast_horizon,
  const double jitter,
  const std::string svd_method
)
{

  // Map the numeric matrices and vectors to eigen objects
  Eigen::Map<Eigen::MatrixXd> X_in_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X_in));
  Eigen::Map<Eigen::VectorXi> delay_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delay));
  Eigen::Map<Eigen::VectorXi> selected_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(selected));
  Eigen::Map<Eigen::VectorXd> l1_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(l1));
  Eigen::Map<Eigen::VectorXd> weights_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(weights));

  // Handle the case where l1, l1_start and or steps is not provided
  if (steps == -2147483647)
  {
    steps = INT_MIN;
  }

  if ((selected_eigen.array() == -2147483647).all())
  {
    selected_eigen.setConstant(INT_MAX);
  }

  if ((l1_eigen.array() == -2147483647).all())
  {
    l1_eigen.setConstant(NAN);
  }

  // Enable/disable parallelisation in Eigen


  // Estimate the sparse DFM
  if (parallel) {
    Eigen::setNbThreads(0);
  }
  else {
    Eigen::setNbThreads(1);
  }
  
  SparseDFM::SDFM<SparseDFM::Structure::SPARSE> results(X_in_eigen, R, order);
  if (svd_method == "fast") {
    results.estimModel<Eigen::BDCSVD>(delay_eigen, selected_eigen, weights_eigen, decorr_errors, crit, l2, l1_eigen, max_iterations, steps, comp_null, spca_conv_crit, fcast_horizon);
  }
  else if (svd_method == "precise") {
    results.estimModel<Eigen::JacobiSVD>(delay_eigen, selected_eigen, weights_eigen, decorr_errors, crit, l2, l1_eigen, max_iterations, steps, comp_null, spca_conv_crit, fcast_horizon);
  }


  // Re-correlate the loadings fit if necessary

  if (decorr_errors)
  {
    Eigen::MatrixXd chol_variable_var_cov = results.inv_chol_variable_var_cov.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(X_in_eigen.cols(), X_in_eigen.cols()));
    results.loading_matrix = (chol_variable_var_cov * results.loading_matrix).eval();
    for (int col = 0; col < results.zero_indeces.cols(); ++col) {
      for (int row = 0; row < results.zero_indeces.rows(); ++row) {
        if (results.zero_indeces(row, col) == 1) {
          results.loading_matrix(row, col) = 0.0;
        }
      }
    }
  }

  Eigen::setNbThreads(0);

  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("Lambda_hat") = Rcpp::wrap(results.loading_matrix),
    Rcpp::Named("Pt") = Rcpp::wrap(results.filter_var_cov),
    Rcpp::Named("F") = Rcpp::wrap(results.factors),
    Rcpp::Named("Wt") = Rcpp::wrap(results.smoother_var_cov),
    Rcpp::Named("C") = Rcpp::wrap(results.inv_chol_variable_var_cov),
    Rcpp::Named("P") = results.max_factor_var_order,
    Rcpp::Named("llt_success_code") = results.llt_success_code);
}

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `twoStepDFM()` instead for a stable and
//' supported interface.
//'
// [[Rcpp::export]]
Rcpp::List runDFMKFS(
  Rcpp::NumericMatrix X_in,
  Rcpp::IntegerVector delay,
  int R,
  int order,
  bool decorr_errors,
  const char* crit,
  double comp_null,
  const bool parallel,
  const unsigned fcast_horizon,
  const double jitter
)
{

  // Map the numeric matrices and vectors to eigen objects
  Eigen::Map<Eigen::MatrixXd> X_in_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X_in));
  Eigen::Map<Eigen::VectorXi> delay_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delay));


  // Enable/disable parallelisation in Eigen


  // Estimate the sparse DFM
  if (parallel) {
    Eigen::setNbThreads(0);
  }
  else {
    Eigen::setNbThreads(1);
  }

  Eigen::VectorXi selected_place_holder = Eigen::VectorXi::Constant(R, X_in_eigen.cols());
  const double l2_placeholder = 0.0;
  Eigen::VectorXd l1_placeholder = INT_MIN * Eigen::VectorXd::Zero(1);
  const int max_iterations_placeholder = 1000;
  int steps_placeholder = INT_MIN;
  const double spca_conv_crit_placeholder = 0.0001;
  Eigen::VectorXd weights_placeholder = Eigen::VectorXd::Ones(X_in_eigen.cols());

  SparseDFM::SDFM<SparseDFM::Structure::DENSE> results(X_in_eigen, R, order);
  // The SVD type does not matter here, as the dense case never calls repeated SVDs.
  results.estimModel<Eigen::JacobiSVD>(delay_eigen, selected_place_holder, weights_placeholder, decorr_errors, crit, l2_placeholder, l1_placeholder, max_iterations_placeholder, steps_placeholder, comp_null, spca_conv_crit_placeholder, fcast_horizon);
  

  // Re-correlate the loadings fit if necessary

  if (decorr_errors)
  {
    Eigen::MatrixXd chol_variable_var_cov = results.inv_chol_variable_var_cov.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(X_in_eigen.cols(), X_in_eigen.cols()));
    results.loading_matrix = (chol_variable_var_cov * results.loading_matrix).eval();
  }

  Eigen::setNbThreads(0);

  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("Lambda_hat") = Rcpp::wrap(results.loading_matrix),
    Rcpp::Named("Pt") = Rcpp::wrap(results.filter_var_cov),
    Rcpp::Named("F") = Rcpp::wrap(results.factors),
    Rcpp::Named("Wt") = Rcpp::wrap(results.smoother_var_cov),
    Rcpp::Named("C") = Rcpp::wrap(results.inv_chol_variable_var_cov),
    Rcpp::Named("P") = results.max_factor_var_order,
    Rcpp::Named("llt_success_code") = results.llt_success_code);
}

/* Simulate an approximate DFM */

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `SimFM()` instead for a stable and
//' supported interface.
//'
// [[Rcpp::export]]
Rcpp::List runStaticFM(
  int T,
  const int& N,
  Rcpp::NumericMatrix S,
  Rcpp::NumericMatrix Lambda,
  Rcpp::NumericVector mu_e,
  Rcpp::NumericMatrix Sigma_e,
  Rcpp::NumericMatrix A,
  int order,
  bool quarterfy,
  bool corr,
  double beta_param,
  double m,
  int seed,
  int R,
  int burn_in,
  bool rescale,
  const bool parallel
)
{

  // Map the numeric matrices and vectors to eigen objects

  Eigen::Map<Eigen::MatrixXd> S_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(S));
  Eigen::Map<Eigen::MatrixXd> Lambda_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Lambda));
  Eigen::Map<Eigen::VectorXd> mu_e_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(mu_e));
  Eigen::Map<Eigen::MatrixXd> Sigma_e_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(Sigma_e));
  Eigen::Map<Eigen::MatrixXd> A_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(A));

  DataGen::FM results;
  std::mt19937 gen(seed);
  if ((burn_in - 1) % 3 == 0) {
    --burn_in;
  }
  else if ((burn_in + 1) % 3 == 0) {
    ++burn_in;
  }

  // Enable/disable parallelisation in Eigen
  if (parallel) {
    Eigen::setNbThreads(0);
  }
  else {
    Eigen::setNbThreads(1);
  }

  DataGen::staticFM(results, T, N, S_eigen, Lambda_eigen, mu_e_eigen, Sigma_e_eigen, A_eigen, gen, order, quarterfy, corr,
    beta_param, m, R, burn_in, rescale);

  Eigen::setNbThreads(0);

  return Rcpp::List::create(Rcpp::Named("F") = Rcpp::wrap(results.F),
    Rcpp::Named("Phi") = Rcpp::wrap(results.Phi),
    Rcpp::Named("Lambda") = Rcpp::wrap(results.Lambda),
    Rcpp::Named("Sigma_xi") = Rcpp::wrap(results.Sigma_e),
    Rcpp::Named("Sigma_epsilon") = Rcpp::wrap(results.Sigma_epsilon),
    Rcpp::Named("Xi") = Rcpp::wrap(results.e),
    Rcpp::Named("X") = Rcpp::wrap(results.X.transpose()),
    Rcpp::Named("frequency") = Rcpp::wrap(results.frequency));

}
