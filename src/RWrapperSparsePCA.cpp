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
Rcpp::List runSPCA(
  Rcpp::NumericMatrix X_in,
  Rcpp::IntegerVector delay,
  Rcpp::IntegerVector selected,
  int R,
  double l2,
  Rcpp::NumericVector l1,
  int max_iterations,
  int steps,
  Rcpp::NumericVector weights,
  double comp_null,
  double spca_conv_crit,
  const bool parallel,
  const std::string svd_method,
  const bool normalise,
  const bool comp_var_expl
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

  int order = 1; // Order does not need to be specified here so its hard-set to 1
  SparseDFM::SDFM<SparseDFM::Structure::SPARSE> results(X_in_eigen, R, order);
  if (svd_method == "fast") { // Use fast BDCSVD for the internal singular value decomposition
    if (comp_var_expl) {
      results.template sparsePrincipalComponents<SparseDFM::SPCAAdditionalComputations::YES, Eigen::BDCSVD>(X_in_eigen.rows() - delay_eigen.maxCoeff(), selected_eigen, l2, l1_eigen, steps, max_iterations, comp_null, spca_conv_crit, normalise, weights_eigen);
    }
    else {
      results.template sparsePrincipalComponents<SparseDFM::SPCAAdditionalComputations::NO, Eigen::BDCSVD>(X_in_eigen.rows() - delay_eigen.maxCoeff(), selected_eigen, l2, l1_eigen, steps, max_iterations, comp_null, spca_conv_crit, normalise, weights_eigen);
    }
  }
  else{
    if (comp_var_expl) { // Use precise JacobiSVD for the internal singular value decomposition
      results.template sparsePrincipalComponents<SparseDFM::SPCAAdditionalComputations::YES, Eigen::JacobiSVD>(X_in_eigen.rows() - delay_eigen.maxCoeff(), selected_eigen, l2, l1_eigen, steps, max_iterations, comp_null, spca_conv_crit, normalise, weights_eigen);
    }
    else {
      results.template sparsePrincipalComponents<SparseDFM::SPCAAdditionalComputations::NO, Eigen::JacobiSVD>(X_in_eigen.rows() - delay_eigen.maxCoeff(), selected_eigen, l2, l1_eigen, steps, max_iterations, comp_null, spca_conv_crit, normalise, weights_eigen);
    }
  }

  Eigen::setNbThreads(0);

  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("Lambda_hat") = Rcpp::wrap(results.loading_matrix),
    Rcpp::Named("F") = Rcpp::wrap(results.factors),
    Rcpp::Named("total_var_expl") = results.total_var_expl,
    Rcpp::Named("pct_var_expl") = results.pct_var_expl
  );
}