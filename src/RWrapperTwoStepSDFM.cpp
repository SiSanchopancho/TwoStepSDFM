/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright \u00A9 2024 Domenic Franjic
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
#include <stdlib.h>
#include <random>
#include <cfloat>
#include <iostream>
#include <Eigen/Eigen>
#include <math.h>
#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "TwoStepSDFM_types.h"

#include <fstream>
#include <string>
#include <chrono>
#include <ctime>

// Internal Incldues
#include "Internals/DataGen.h" // Data generation
#include "Internals/SparseDFM.h" // Wrapper for the Sparse DFM estimation

/* Function to run the two-step SDFM estimation procedure*/
/*Source:
 -  Franjic, Domenic and Schweikert, Karsten, Nowcasting Macroeconomic Variables with a Sparse Mixed Frequency Dynamic Factor Model (February 21, 2024). Available at SSRN: https://ssrn.com/abstract=4733872 or http://dx.doi.org/10.2139/ssrn.4733872
 */

//' @description
//' This function is for internal use only and may change in future releases
//' without notice. Users should use `SimFM()` instead for a stable and
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
  double comp_null,
  bool check_rank,
  double conv_crit,
  double conv_threshold,
  bool log,
  double KFS_conv_crit,
  const bool parallel,
  const unsigned fcast_horizon
)
{

  // Initialise the result object
  Filtering::KFS_fit results;

  // Map the numeric matrices and vectors to eigen objects
  Eigen::Map<Eigen::MatrixXd> X_in_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X_in));
  Eigen::Map<Eigen::VectorXi> delay_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(delay));
  Eigen::Map<Eigen::VectorXi> selected_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXi>>(selected));
  Eigen::Map<Eigen::VectorXd> l1_eigen(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(l1));


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

  SparseDFM::SDFMKFS(results, X_in_eigen, delay_eigen, selected_eigen, R, order, decorr_errors, crit,
    l2, l1_eigen, max_iterations, steps, comp_null, check_rank, conv_crit, conv_threshold, log,
    KFS_conv_crit, fcast_horizon);

  // Re-correlate the loadings fit if necessary

  if (decorr_errors)
  {
    Eigen::MatrixXd C_inv = results.C.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(X_in_eigen.cols(), X_in_eigen.cols()));
    results.Lambda_hat = (C_inv * results.Lambda_hat).eval();
    for (int col = 0; col < results.Zero_Indeces.cols(); ++col) {
      for (int row = 0; row < results.Zero_Indeces.rows(); ++row) {
        if (results.Zero_Indeces(row, col) == 1) {
          results.Lambda_hat(row, col) = 0.0;
        }
      }
    }
  }

  Eigen::setNbThreads(0);

  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("Lambda_hat") = Rcpp::wrap(results.Lambda_hat),
    Rcpp::Named("Pt") = Rcpp::wrap(results.Pt),
    Rcpp::Named("F") = Rcpp::wrap(results.F),
    Rcpp::Named("Wt") = Rcpp::wrap(results.Wt),
    Rcpp::Named("C") = Rcpp::wrap(results.C),
    Rcpp::Named("P") = results.order);

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