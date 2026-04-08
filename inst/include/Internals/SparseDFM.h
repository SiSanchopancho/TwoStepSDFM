/* SPDX-License-Identifier: GPL-3.0-or-later */
/*
 * Copyright (cholesky_factor_variable_var_cov) 2024 Domenic Franjic
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

#pragma once

#ifndef SPARSE_DYNAMIC_FACTOR_MODEL
#define SPARSE_DYNAMIC_FACTOR_MODEL

#define _USE_MATH_DEFINES

 // Including external libraries
#if !defined(_MSC_VER)
#define EIGEN_NO_DEBUG
#include <RcppEigen.h>
#else
#include <Eigen>
#include "Developer.h"
#endif
#include <stdlib.h>
#include <string>
#include <math.h>

// Including internal libraries
#include "Orders.h"
#include "ElNetSolve.h"

namespace SparseDFM {

  enum Structure {
    SPARSE = 1, DENSE = 2
  };

  enum SPCAAdditionalComputations {
    YES = 1, NO = 2
  };

  template <Structure structure_type>
  class SDFM {
  public:
    Eigen::MatrixXd data;
    Eigen::MatrixXd loading_matrix;
    Eigen::MatrixXi zero_indeces;
    Eigen::MatrixXd filter_var_cov;
    Eigen::MatrixXd factors;
    Eigen::MatrixXd smoother_var_cov;
    Eigen::MatrixXd inv_chol_variable_var_cov;
    int no_of_factors;
    int no_of_obs;
    int no_of_vars;
    int max_factor_var_order;
    bool conv;
    double total_var_expl;
    Eigen::VectorXd pct_var_expl;
    int llt_success_code;

    // Default Constructor
    SDFM() : no_of_factors(0), no_of_obs(0), no_of_vars(0), max_factor_var_order(0), conv(true), total_var_expl(0.0), llt_success_code(0){
      this->data = Eigen::MatrixXd::Zero(0, 0);
      this->loading_matrix = Eigen::MatrixXd::Zero(0, 0);
      this->zero_indeces = Eigen::MatrixXi::Zero(0, 0);
      this->filter_var_cov = Eigen::MatrixXd::Zero(0, 0);
      this->factors = Eigen::MatrixXd::Zero(0, 0);
      this->smoother_var_cov = Eigen::MatrixXd::Zero(0, 0);
      this->inv_chol_variable_var_cov = Eigen::MatrixXd::Zero(0, 0);
      this->pct_var_expl = Eigen::VectorXd::Zero(0);
    }
    SDFM(const Eigen::MatrixXd& data, int no_of_factors, int max_factor_var_order) : no_of_factors(no_of_factors), max_factor_var_order(max_factor_var_order), conv(false), total_var_expl(0.0), llt_success_code(0) {
      this->data = data;
      this->no_of_vars = data.cols();
      this->no_of_obs = data.rows();
      this->loading_matrix = Eigen::MatrixXd::Zero(this->no_of_vars, this->no_of_factors);
      this->zero_indeces = Eigen::MatrixXi::Zero(this->no_of_vars, this->no_of_factors);
      this->filter_var_cov = Eigen::MatrixXd::Zero(this->no_of_factors, this->no_of_factors);
      this->factors = Eigen::MatrixXd::Zero(this->no_of_factors, this->no_of_obs);
      this->smoother_var_cov = Eigen::MatrixXd::Zero(this->no_of_factors, this->no_of_factors);
      this->inv_chol_variable_var_cov = Eigen::MatrixXd::Zero(this->no_of_vars, this->no_of_vars);
      this->pct_var_expl = Eigen::VectorXd::Zero(0); // This is only initilaised in sparsePrincipalComponentAnalysis if its computation is explicitly required
    }
    SDFM(const SDFM& other) = default;
    SDFM(SDFM&& other) noexcept = default;
    SDFM& operator=(const SDFM& other) = default;
    SDFM& operator=(SDFM&& other) noexcept = default;
    ~SDFM() = default;

    /* Main estimation funciton for the two-step sparse dnymamic factor model */

    template <template <typename> class SVDType>
    void estimModel(
      const Eigen::VectorXi& delay,
      const Eigen::VectorXi& selected,
      const Eigen::VectorXd& weights,
      const bool& decorr_errors = 0,
      const char* crit = "BIC",
      const double& ridge = 0.4,
      Eigen::VectorXd lasso_penalties = INT_MIN * Eigen::VectorXd::Zero(1),
      const int& max_iterations = 1000,
      int steps = INT_MIN,
      const double& comp_null = 10e-15,
      const double& spca_conv_crit = 0.0001,
      const unsigned& fcast_horizon = 0,
      const double& jitter = 1e-8
    ) {

      /* Compute the SPCA estimates */

      Eigen::MatrixXd effective_data = this->data.transpose();
      if (structure_type == Structure::SPARSE) {
        const bool normalise = false;
        sparsePrincipalComponents<SPCAAdditionalComputations::NO, SVDType>(this->no_of_obs - delay.maxCoeff(), selected, ridge, lasso_penalties, steps, max_iterations, comp_null, spca_conv_crit, normalise, weights);
      }
      else if (structure_type == Structure::DENSE) {
        principalComponents(this->no_of_obs - delay.maxCoeff());
      }
      int factor_var_order = VARorder<BIC>(this->factors, this->max_factor_var_order, comp_null);
      int ind_of_first_companion_obs = factor_var_order - 1;


      /* Initialise the kalman filter */

      // Build companion form matrices and cut data to size
      int no_of_states = this->no_of_factors * factor_var_order;
      Eigen::MatrixXd comp_form_factors = Eigen::MatrixXd::Zero(this->no_of_factors * factor_var_order, this->no_of_obs - delay.maxCoeff());
      Eigen::MatrixXd filter_var_cov_initial = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
      filter_var_cov_initial.diagonal().setConstant(1000);
      Eigen::MatrixXd comp_form_loading_matrix = Eigen::MatrixXd::Zero(this->no_of_vars, no_of_states);
      Eigen::MatrixXd state_equation_var_cov = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
      Eigen::MatrixXd meas_equation_var_cov = Eigen::MatrixXd::Zero(this->no_of_vars, this->no_of_vars);
      Eigen::MatrixXd comp_form_factor_var_coeff = Eigen::MatrixXd::Zero(no_of_states, no_of_states);

      for (int curr_lag = 0; curr_lag < factor_var_order; ++curr_lag) {
        comp_form_factors.block(this->no_of_factors * curr_lag, curr_lag, this->no_of_factors, this->no_of_obs - delay.maxCoeff() - curr_lag) = this->factors.block(0, 0, no_of_factors, no_of_obs - delay.maxCoeff() - curr_lag);
      }

      // Estimate the missing state and measurement equation parameters
      int index_of_full_sample = no_of_obs - delay.maxCoeff() - 1;
      {

        // State equation parameters
        Eigen::MatrixXd current_comp_form_factors = comp_form_factors(Eigen::all, Eigen::seq(ind_of_first_companion_obs + 1, Eigen::last)).transpose();
        Eigen::MatrixXd lagged_comp_form_factors = comp_form_factors(Eigen::all, Eigen::seq(ind_of_first_companion_obs, Eigen::last - 1)).transpose();
        Eigen::MatrixXd ident_no_of_states = Eigen::MatrixXd::Identity(no_of_states, no_of_states);
        comp_form_factor_var_coeff = ((lagged_comp_form_factors.transpose() * lagged_comp_form_factors).llt().solve(ident_no_of_states) * lagged_comp_form_factors.transpose() * current_comp_form_factors).transpose();
        comp_form_factor_var_coeff.bottomLeftCorner(this->no_of_factors * (ind_of_first_companion_obs), this->no_of_factors * (ind_of_first_companion_obs)).setIdentity();
        comp_form_factor_var_coeff.bottomRightCorner(this->no_of_factors * (ind_of_first_companion_obs), this->no_of_factors).setZero();
        Eigen::MatrixXd factor_var_coeff_transposed = comp_form_factor_var_coeff.transpose();
        Eigen::MatrixXd state_equation_residuals = current_comp_form_factors - lagged_comp_form_factors * factor_var_coeff_transposed;
        //state_equation_var_cov.diagonal() = (1. / double(this->no_of_obs - factor_var_order - 1)) * state_equation_residuals.array().square().colwise().sum();
        state_equation_var_cov.topLeftCorner(this->no_of_factors, this->no_of_factors) = (1. / double(this->no_of_obs - factor_var_order - 1)) * (state_equation_residuals * state_equation_residuals.transpose()).topLeftCorner(this->no_of_factors, this->no_of_factors);

        // Measurement equation parameters
        comp_form_loading_matrix.topLeftCorner(this->no_of_vars, this->no_of_factors) = this->loading_matrix;
        Eigen::MatrixXd meas_equation_residuals = effective_data(Eigen::all, Eigen::seq(0, index_of_full_sample)) - comp_form_loading_matrix * comp_form_factors(Eigen::all, Eigen::seq(0, index_of_full_sample));
        if (decorr_errors) {
          meas_equation_var_cov = meas_equation_residuals * meas_equation_residuals.transpose();
        }
        else {
          meas_equation_var_cov.diagonal() = meas_equation_residuals.array().square().colwise().sum();
        }
        meas_equation_var_cov *= (1. / double(this->no_of_obs - delay.maxCoeff()));
        meas_equation_var_cov.diagonal().array() += jitter; // jitter the matrix for numerical stability
      }
      this->factors.conservativeResize(no_of_states, (this->no_of_obs + fcast_horizon + 1));
      this->factors.setZero();
      this->filter_var_cov.conservativeResize(no_of_states, no_of_states * (this->no_of_obs + fcast_horizon + 1));
      this->filter_var_cov.setZero();
      this->smoother_var_cov.conservativeResize(no_of_states, no_of_states * this->no_of_obs);
      this->smoother_var_cov.setZero();
      this->zero_indeces.conservativeResize(this->loading_matrix.rows(), this->loading_matrix.cols());
      this->zero_indeces.setZero();

      if (decorr_errors) { // Decorrelate the data and filter dynamics for the implementation of the univariate Kalman filter and smoother
        this->zero_indeces = this->loading_matrix.unaryExpr([comp_null](double a) { return (a == 0) ? 1 : 0; });
        Eigen::LLT<Eigen::MatrixXd> meas_llt(meas_equation_var_cov);
        if (meas_llt.info() != Eigen::Success) {
#if !defined(_MSC_VER)
          this->llt_success_code = -1;
#else
          std::cerr << "\nWARNING: Decorrelation failed. Using robust cholesky decomposition.\n";
#endif
          Eigen::LDLT<Eigen::MatrixXd> meas_ldlt(meas_equation_var_cov);
          if (meas_ldlt.info() != Eigen::Success) {
#if !defined(_MSC_VER)
            this->llt_success_code = -2;
#else
            std::cerr << "\nWARNING: Robust decorrelation failed. Ignoring measurement error correlation structure.\n";
#endif
            this->inv_chol_variable_var_cov = Eigen::MatrixXd::Identity(this->no_of_vars, this->no_of_vars);
          }
          Eigen::MatrixXd lower_factor = meas_ldlt.matrixL();
          Eigen::VectorXd diagonal_matrix_diag = meas_ldlt.vectorD();
          diagonal_matrix_diag.diagonal() = diagonal_matrix_diag.diagonal().unaryExpr([comp_null](double x) { return std::max(x, comp_null); }).eval();
          Eigen::DiagonalMatrix<double, Eigen::Dynamic> sqrt_diagonal_matrix_diag = meas_ldlt.vectorD().cwiseSqrt().asDiagonal();
          Eigen::Transpositions<Eigen::Dynamic> permutation_matrix = meas_ldlt.transpositionsP();
          Eigen::MatrixXd chol_variable_var_cov = lower_factor * sqrt_diagonal_matrix_diag;
          this->inv_chol_variable_var_cov = chol_variable_var_cov.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(this->no_of_vars, this->no_of_vars));;
        }
        else {
          this->inv_chol_variable_var_cov = meas_llt.matrixL().solve(Eigen::MatrixXd::Identity(this->no_of_vars, this->no_of_vars));
        }

        if (0 < delay.maxCoeff()) { // Decorrelate row by row in the case of missing data
          Eigen::MatrixXd Temp = effective_data;
          for (int curr_obs = 0; curr_obs < this->no_of_obs; ++curr_obs) {
            Eigen::ArrayXd missing_ind = (!((Eigen::VectorXi::Constant(this->no_of_vars, 1, this->no_of_obs) - delay).array() <= curr_obs)).template cast<double>();
            for (int curr_var = 0; curr_var < this->no_of_vars; ++curr_var)
            {
              if (this->no_of_obs - delay(curr_var) <= curr_obs) {
                Temp(curr_var, curr_obs) = DBL_MAX;
              }
              else {
                Eigen::VectorXd Cn = (this->inv_chol_variable_var_cov.row(curr_var).transpose().array() * missing_ind).matrix();
                Eigen::VectorXd Xt = effective_data.col(curr_obs).array().isNaN().select(0.0, effective_data.col(curr_obs)).matrix();
                Temp(curr_var, curr_obs) = Cn.dot(Xt);
              }
            }
          }
          effective_data = Temp;
        }
        else
        {
          effective_data = this->inv_chol_variable_var_cov * effective_data;
        }
        comp_form_loading_matrix.topLeftCorner(this->no_of_vars, this->no_of_factors) = (this->inv_chol_variable_var_cov * comp_form_loading_matrix.topLeftCorner(this->no_of_vars, this->no_of_factors)).eval();
        meas_equation_var_cov.setZero();
        meas_equation_var_cov.diagonal() = (1. / double(this->no_of_obs - delay.maxCoeff()))
          * (effective_data(Eigen::all, Eigen::seq(0, index_of_full_sample))
            - comp_form_loading_matrix
            * comp_form_factors(Eigen::all, Eigen::seq(0, index_of_full_sample))
            ).array().square().rowwise().sum();

      }
      else { // Store identity as inverse cholesky factor in the case of no decorrelation
        this->inv_chol_variable_var_cov = Eigen::MatrixXd::Identity(this->no_of_vars, this->no_of_vars);
      }

      /* State filtering and smoothing */

      Eigen::MatrixXd inv_meas_equation_var_cov = Eigen::MatrixXd::Zero(this->no_of_vars, this->no_of_obs);
      Eigen::MatrixXd kalman_gain = Eigen::MatrixXd::Zero(no_of_states, this->no_of_vars * this->no_of_obs);
      Eigen::MatrixXd filter_error = Eigen::MatrixXd::Zero(this->no_of_vars, this->no_of_obs);

      int kalman_filter_starting_index = ind_of_first_companion_obs;
      univariateRepOfMultivariateKalmanFilter(comp_form_factors, inv_meas_equation_var_cov, kalman_gain, filter_error, kalman_filter_starting_index, this->no_of_obs, no_of_states, effective_data, comp_form_loading_matrix, comp_form_factor_var_coeff, meas_equation_var_cov, state_equation_var_cov, filter_var_cov_initial, delay, fcast_horizon);
      univariateRepOfMultivariateKalmanSmoother(kalman_gain, inv_meas_equation_var_cov, filter_error, this->no_of_obs, no_of_states, effective_data, comp_form_loading_matrix, comp_form_factor_var_coeff, meas_equation_var_cov, state_equation_var_cov, delay);
      this->loading_matrix = comp_form_loading_matrix;
      this->max_factor_var_order = factor_var_order;

      return;
    }

    /* Loading matrix estimation routines */

    void principalComponents(const int effective_time) {

      // PCA using the Eigen decomposition
      Eigen::RowVectorXd mean = this->data(Eigen::seq(0, effective_time - 1), Eigen::all).array().colwise().mean();
      Eigen::MatrixXd centered_data = this->data(Eigen::seq(0, effective_time - 1), Eigen::all) - mean.replicate(effective_time, 1);
      Eigen::MatrixXd data_var_cov = centered_data.transpose() * centered_data * 1.0 / (effective_time - 1.0);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_deco(data_var_cov);
      this->loading_matrix = eigen_deco.eigenvectors().rightCols(this->no_of_factors); // Use the most right columns, as Eigen stores eigenvaluers from smallest to largest
      for (int factor = 0; factor < this->no_of_factors; ++factor)
      {
        this->loading_matrix.col(factor).normalize();
      }
      this->factors = (this->data(Eigen::seq(0, effective_time - 1), Eigen::all) * this->loading_matrix).transpose();

      return;
    }

    /* Sparse Principal Components Analysis */
    /* Source:
        - Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. Journal of computational and graphical statistics, 15(2), 265-286. (https://doi.org/10.1198/106186006X113430)
        - Zou, H., Hastie, T., & Zou, M. H. (2016). Package ?elasticnet?. (https://cran.r-project.org/web/packages/elasticnet/index.html)
        - Zou, Q., & Zhang, P. (2024, December). On General Weighted Adaptive Sparse Principal Component Analysis. In Proceedings of the 2024 4th International Conference on Computational Modeling, Simulation and Data Analysis (pp. 335-340).
    */
    template <SPCAAdditionalComputations compute_add_stuff, template <typename> class SVDType> // The type of SVD is templated to give the user the choice between using Eigen's save JacobiSVD agains fast BDCSVD
    void sparsePrincipalComponents(
      const int effective_time,
      const Eigen::VectorXi& selected,
      const double& ridge,
      Eigen::VectorXd lasso_penalties,
      int steps,
      const int max_iterations,
      const double& comp_null, 
      const double& spca_conv_crit, 
      const bool& normalise,
      const Eigen::VectorXd& weights
    ) {

      // Compute initial data SVD and store its matrices
      SVDType<Eigen::MatrixXd> svd;
      svd.compute(this->data(Eigen::seq(0, effective_time - 1), Eigen::all), Eigen::ComputeThinU | Eigen::ComputeThinV);
      if (compute_add_stuff == SPCAAdditionalComputations::YES) {
        this->total_var_expl = svd.singularValues().array().square().sum();
      }

      const Eigen::MatrixXd effective_data = data(Eigen::seq(0, effective_time - 1), Eigen::all);
      Eigen::MatrixXd gram = effective_data.transpose() * effective_data;
      Eigen::VectorXd artificial_target = Eigen::VectorXd::Zero(effective_time);
      Eigen::MatrixXd data_lambda_gram = Eigen::MatrixXd::Zero(this->no_of_vars, this->no_of_vars);
      Eigen::MatrixXd data_svd_u_mat = svd.matrixU();
      Eigen::MatrixXd data_svd_v_mat = svd.matrixV();
      Eigen::MatrixXd sing_val = svd.singularValues().asDiagonal();

      /* Start initial SPCA block */

      Eigen::MatrixXd dual_matrix = data_svd_v_mat.leftCols(this->no_of_factors);
      Eigen::MatrixXd past_dual_matrix = dual_matrix;
      const bool pure_ridge_case = (steps == INT_MIN) && (selected.array() >= this->no_of_vars).all() && lasso_penalties.array().isNaN().any();
      if (pure_ridge_case) { // This represents the pure ridge case. No LARS needed
        Eigen::LLT<Eigen::MatrixXd> jittered_gram_llt(gram + Eigen::MatrixXd::Identity(this->no_of_vars, this->no_of_vars) * ridge);
        for (int factor = 0; factor < this->no_of_factors; ++factor) {
          artificial_target = effective_data * dual_matrix.col(factor);
          this->loading_matrix.col(factor) = jittered_gram_llt.solve(effective_data.transpose() * artificial_target);
        }
      }
      else { // Use LARS-EN to solve the LS problem with size side-constraints
        for (int factor = 0; factor < this->no_of_factors; ++factor) {
          artificial_target = effective_data * dual_matrix.col(factor);
          this->loading_matrix.col(factor) = LARS<false>(artificial_target, effective_data, weights, ridge, lasso_penalties(factor), selected(factor), steps, comp_null);
        }
      }
      
      /* End initial SPCA block */

      /* Start SPCA refinement loop */

      for (int i = 0; i < max_iterations; ++i) {

        data_lambda_gram = gram * this->loading_matrix;
        svd.compute(data_lambda_gram, Eigen::ComputeThinU | Eigen::ComputeThinV);
        dual_matrix = svd.matrixU() * svd.matrixV().transpose();

        if ((past_dual_matrix - dual_matrix).squaredNorm() <= spca_conv_crit) { // Convergence check
          break;
        }
        else {
          past_dual_matrix = dual_matrix;
        }

        if (steps == INT_MIN && (selected.array() >= this->no_of_vars).all() && lasso_penalties.array().isNaN().any()) { // This represents the pure ridge case. No LARS needed
          Eigen::LLT<Eigen::MatrixXd> jittered_gram_llt(gram + Eigen::MatrixXd::Identity(this->no_of_vars, this->no_of_vars) * ridge);
          for (int factor = 0; factor < this->no_of_factors; ++factor) {
            artificial_target = effective_data * dual_matrix.col(factor);
            this->loading_matrix.col(factor) = jittered_gram_llt.solve(effective_data.transpose() * artificial_target);
          }
        }
        else { // Use LARS-EN to solve the LS problem with size side-constraints
          for (int factor = 0; factor < this->no_of_factors; ++factor) {
            artificial_target = effective_data * dual_matrix.col(factor);
            this->loading_matrix.col(factor) = LARS<false>(artificial_target, effective_data, weights, ridge, lasso_penalties(factor), selected(factor), steps, comp_null);
          }
        }
      }

      /* End SPCA refinement loop */

      // Calculate the factors
      this->factors = (effective_data * this->loading_matrix).transpose();

      if (compute_add_stuff == SPCAAdditionalComputations::YES) { // Calculate variance explained if additional information is asked for
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR;
        Eigen::MatrixXd R = QR.compute(this->factors.transpose()).matrixQR().template triangularView<Eigen::Upper>();
        this->pct_var_expl= 1 / this->total_var_expl * (R.diagonal().array().square()).matrix();
      }

      if (normalise) {
        for (int factor = 0; factor < this->no_of_factors; ++factor)
        {
          this->loading_matrix.col(factor).normalize();
        }
      }

      return;

    }

    /* Filter routines */

    void univariateRepOfMultivariateKalmanFilter(
      Eigen::MatrixXd& comp_form_factors,
      Eigen::MatrixXd& inv_meas_equation_var_cov,
      Eigen::MatrixXd& kalman_gain,
      Eigen::MatrixXd& filter_error,
      const int& starting_index,
      const int& no_of_obs,
      const int& no_of_states,
      const Eigen::MatrixXd& effective_data,
      const Eigen::MatrixXd& comp_form_loading_matrix,
      const Eigen::MatrixXd& comp_form_factor_var_coeff,
      const Eigen::MatrixXd& meas_equation_var_cov,
      const Eigen::MatrixXd& state_equation_var_cov,
      const Eigen::MatrixXd& filter_var_cov_initial,
      const Eigen::VectorXi& delay,
      const int& forecast_horizon
    )
    {

      Eigen::MatrixXd state_initial = comp_form_factors.rowwise().mean();
      double current_variable_var_inv = 0.0;
      Eigen::VectorXd current_state_estimate = Eigen::VectorXd::Zero(no_of_states);
      Eigen::MatrixXd current_uncertainty = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
      Eigen::MatrixXd comp_form_loading_matrix_transopose = comp_form_loading_matrix.transpose();
      this->factors.leftCols(starting_index + 1) = state_initial.replicate(1, starting_index + 1);
      this->filter_var_cov.topLeftCorner(no_of_states, (starting_index + 1) * no_of_states) = filter_var_cov_initial.replicate(1, starting_index + 1);

      /* Start loop over observarions */

      for (int curr_obs = 0; curr_obs < no_of_obs; ++curr_obs) {
        current_state_estimate = this->factors.col(curr_obs);
        current_uncertainty = this->filter_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states);

        /* Start loop over variables */

        for (int curr_var = 0; curr_var < this->no_of_vars; ++curr_var) {
          if (no_of_obs - delay(curr_var) <= curr_obs) { // Skip missing observations
            continue;
          }

          /* Update */

          current_variable_var_inv = comp_form_loading_matrix.row(curr_var) * current_uncertainty * comp_form_loading_matrix_transopose.col(curr_var) + meas_equation_var_cov.diagonal()(curr_var);
          filter_error(curr_var, curr_obs) = effective_data(curr_var, curr_obs) - comp_form_loading_matrix.row(curr_var) * current_state_estimate;
          inv_meas_equation_var_cov(curr_var, curr_obs) = 1.0 / current_variable_var_inv;
          int curr_kalman_gain_ind = curr_obs * this->no_of_vars + curr_var;
          kalman_gain.col(curr_kalman_gain_ind) = current_uncertainty * comp_form_loading_matrix_transopose.col(curr_var);
          current_state_estimate += (inv_meas_equation_var_cov(curr_var, curr_obs)) * kalman_gain.col(curr_kalman_gain_ind) * filter_error(curr_var, curr_obs);
          current_uncertainty -= kalman_gain.col(curr_kalman_gain_ind) * (inv_meas_equation_var_cov(curr_var, curr_obs)) * kalman_gain.col(curr_kalman_gain_ind).transpose();

        }

        /* End loop over variables */

        /* Store and propagade */

        this->factors.col(curr_obs) = current_state_estimate;
        this->filter_var_cov.block(0, (curr_obs)*no_of_states, no_of_states, no_of_states) = current_uncertainty;
        this->factors.col(curr_obs + 1) = comp_form_factor_var_coeff * current_state_estimate;
        this->filter_var_cov.block(0, (curr_obs + 1) * no_of_states, no_of_states, no_of_states) = comp_form_factor_var_coeff * current_uncertainty * comp_form_factor_var_coeff.transpose() + state_equation_var_cov;
      }

      /* End loop over observarions */

      /* Forecast */

      if (0 < forecast_horizon) {
        for (int curr_obs = no_of_obs; curr_obs < no_of_obs + forecast_horizon; ++curr_obs) {
          this->factors.col(curr_obs + 1) = comp_form_factor_var_coeff * this->factors.col(curr_obs);
          this->filter_var_cov.block(0, (curr_obs + 1) * no_of_states, no_of_states, no_of_states) = comp_form_factor_var_coeff * this->filter_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states) * comp_form_factor_var_coeff.transpose() + state_equation_var_cov;
        }
      }
      return;
    }

    /* Univariate representation of the Multivariate Kalman Smoother */
    /* Sources:
    Koopman, S. J., & Durbin, J. (2000). Fast filtering and smoothing for multivariate state space models. Journal of time series analysis, 21(3), 281-296. (https://doi.org/10.1111/1467-9892.00186)
    */
    void univariateRepOfMultivariateKalmanSmoother(
      Eigen::MatrixXd& kalman_gain,
      Eigen::MatrixXd& inv_meas_equation_var_cov,
      Eigen::MatrixXd& filter_error,
      const int& no_of_obs,
      const int& no_of_states,
      const Eigen::MatrixXd& effective_data,
      const Eigen::MatrixXd& comp_form_loading_matrix,
      const Eigen::MatrixXd& comp_form_factor_var_coeff,
      const Eigen::MatrixXd& meas_equation_var_cov,
      const Eigen::MatrixXd& state_equation_var_cov,
      const Eigen::VectorXi& delay
    )
    {

      Eigen::VectorXd current_smoothed_state = Eigen::VectorXd::Zero(no_of_states);
      Eigen::MatrixXd weighted_kalman_gain = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
      Eigen::MatrixXd current_smoothed_uncert = Eigen::MatrixXd::Zero(no_of_states, no_of_states);
      Eigen::MatrixXd comp_form_loading_matrix_transopose = comp_form_loading_matrix.transpose();
      Eigen::MatrixXd ident_no_of_states = Eigen::MatrixXd::Identity(no_of_states, no_of_states);

      /* Start loop over observations */

      for (int curr_obs = no_of_obs - 1; curr_obs >= 0; --curr_obs) {

        /* Start loop over variables */

        for (int curr_var = this->no_of_vars - 1; curr_var >= 0; --curr_var) {

          if (no_of_obs - delay(curr_var) <= curr_obs) { // Skip missing observations
            continue;
          }

          /* update */

          weighted_kalman_gain = ident_no_of_states - kalman_gain.col(curr_obs * this->no_of_vars + curr_var) * comp_form_loading_matrix.row(curr_var) * inv_meas_equation_var_cov(curr_var, curr_obs);
          current_smoothed_state = (comp_form_loading_matrix_transopose.col(curr_var) * inv_meas_equation_var_cov(curr_var, curr_obs) * filter_error(curr_var, curr_obs) + weighted_kalman_gain.transpose() * current_smoothed_state).eval();
          current_smoothed_uncert = (comp_form_loading_matrix_transopose.col(curr_var) * inv_meas_equation_var_cov(curr_var, curr_obs) * comp_form_loading_matrix.row(curr_var) + weighted_kalman_gain.transpose() * current_smoothed_uncert * weighted_kalman_gain).eval();

        }

        /* End loop over variables */

        /* Smooth and propagade */
        this->factors.col(curr_obs) += this->filter_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states) * current_smoothed_state;
        this->smoother_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states) = this->filter_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states) - this->filter_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states) * current_smoothed_uncert * this->filter_var_cov.block(0, curr_obs * no_of_states, no_of_states, no_of_states);
        current_smoothed_state = (comp_form_factor_var_coeff.transpose() * current_smoothed_state).eval();
        current_smoothed_uncert = (comp_form_factor_var_coeff.transpose() * current_smoothed_uncert * comp_form_factor_var_coeff).eval();

      }

      /* End loop over observations */

      return;
    }

  };

};
#endif /* defined(SPARSE_DYNAMIC_FACTOR_MODEL) */

