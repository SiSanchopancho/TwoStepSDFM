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

#if !defined(_MSC_VER)
 // Source: Onatski, A. (2009). Supplement to "Testing hypotheses about the number of factors in large factor models". In Econometrica, 77(5). The Econometric Society. https://www.econometricsociety.org/publications/econometrica/issue-supplemental-materials/2009/09
  //' @description
   //' This function is for internal use only and may change in future releases
   //' without notice. 
   //'
   // [[Rcpp::export]]
Rcpp::List runNoOfFactorsTest(
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
    Rcpp::Named("confidence_threshold") = Rcpp::wrap(confidence_threshold),
    Rcpp::Named("eigen_values") = Rcpp::wrap(eigen_values)
  );

}

// Source: Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. Econometrica, 70(1), 191-221.
 //' @description
  //' This function is for internal use only and may change in future releases
  //' without notice. 
  //'
  // [[Rcpp::export]]
Rcpp::List runNoOfFactorsInfoCrit(
  Rcpp::NumericMatrix data_matrix_in,
  const int max_no_factors
)
{

  // PCA using the Eigen decomposition
  Eigen::Map<Eigen::MatrixXd> data(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(data_matrix_in));
  Eigen::VectorXd mean = data.array().rowwise().mean();
  Eigen::MatrixXd centered_data = data - mean.replicate(1, data.cols());
  Eigen::MatrixXd data_var_cov = centered_data * centered_data.transpose() * (1.0 / (data.cols() - 1.0));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_deco(data_var_cov);
  Eigen::MatrixXd loading_matrix = eigen_deco.eigenvectors().rowwise().reverse().leftCols(max_no_factors);
  const int no_of_vars = static_cast<double>(data.rows());
  const int no_of_obs = static_cast<double>(data.cols());
  for (int factor = 0; factor < max_no_factors; ++factor)
  {
    loading_matrix.col(factor).normalize();
  }
  Eigen::MatrixXd factors = loading_matrix.transpose() * data;

  Eigen::MatrixXd information_crit = Eigen::MatrixXd(max_no_factors, 3);
  for (int curr_fact = 0; curr_fact < max_no_factors; ++curr_fact) {
    Eigen::MatrixXd residuals = data - loading_matrix.leftCols(curr_fact + 1) * factors.topRows(curr_fact + 1);
    const double dim_penalty = std::min(no_of_vars, no_of_obs);
    const double common_penalty = static_cast<double>(curr_fact) * ((no_of_vars + no_of_obs) / (no_of_vars * no_of_obs));
    const double log_rss = std::log(residuals.array().square().mean());
    information_crit(curr_fact, 0) = log_rss + common_penalty * log((no_of_vars * no_of_obs) / (no_of_vars + no_of_obs));
    information_crit(curr_fact, 1) = log_rss + common_penalty * log(dim_penalty);
    information_crit(curr_fact, 2) = log_rss + static_cast<double>(curr_fact) * (std::log(dim_penalty) / dim_penalty);
  }
  // Convert the results back to Rcpp types and return
  return Rcpp::List::create(Rcpp::Named("information_crit") = Rcpp::wrap(information_crit)
  );

}
#endif