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

#include "../inst/include/Internals/Orders.h"

// Source: Onatski, A. (2009). Supplement to "Testing hypotheses about the number of factors in large factor models". In Econometrica, 77(5). The Econometric Society. https://www.econometricsociety.org/publications/econometrica/issue-supplemental-materials/2009/09
 //' @description
  //' This function is for internal use only and may change in future releases
  //' without notice. 
  //'
  // [[Rcpp::export]]
Rcpp::List runNoOfFactors(
  Rcpp::NumericMatrix data_matrix_in,
  Rcpp::NumericMatrix test_values,
  const int min_no_factors,
  const int max_no_factors,
  const double confidence_threshold
)
{

  Eigen::Map<Eigen::MatrixXd> data(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(data_matrix_in));
  Eigen::Map<Eigen::MatrixXd> test_values_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(test_values));
  int no_of_factors = min_no_factors - 1, no_of_obs = data.cols(), cutoff = std::floor(double(no_of_obs) / 2.);
  double p_value = DBL_MIN, test_statistic = DBL_MIN;
  std::complex<double> i(0, 1);

  Eigen::MatrixXcd data_complex = data(Eigen::all, Eigen::seq(0, cutoff - 1)) + i * data(Eigen::all, Eigen::seq(cutoff, 2 * cutoff - 1));

  // Eigen decompositions
  Eigen::MatrixXcd gram = data_complex * data_complex.adjoint();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(gram);
  Eigen::VectorXd eigen_values = eig.eigenvalues().reverse();

  /* Start of the test loop */

  while (p_value < confidence_threshold && no_of_factors < max_no_factors - 1)
  {
    
    ++no_of_factors;

    // Calculate the test statistic
    test_statistic = ((eigen_values(Eigen::seq(no_of_factors, max_no_factors - 1)) - eigen_values(Eigen::seq(no_of_factors + 1, max_no_factors))).array() / (eigen_values(Eigen::seq(no_of_factors + 1, max_no_factors)) - eigen_values(Eigen::seq(no_of_factors + 2, max_no_factors + 1))).array()).maxCoeff();
    p_value = static_cast<double>((test_values_eigen(Eigen::all, max_no_factors - no_of_factors - 1).array() > test_statistic).count()) / 1000.0;

  }

  /* End of the test loop */

  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("no_of_factors") = Rcpp::wrap(no_of_factors),
    Rcpp::Named("p_value") = Rcpp::wrap(p_value),
    Rcpp::Named("test_statistic") = Rcpp::wrap(test_statistic),
    Rcpp::Named("confidence_threshold") = Rcpp::wrap(confidence_threshold)
  );

}