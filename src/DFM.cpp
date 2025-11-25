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


#include "Internals/SparseDFM.h"

 /* Two-Step Estimator for a LInear Gaussian Sparse Latent (Static) Dynamic Factor Model */
void SparseDFM::DFMKFS(
	Filtering::KFS_fit& results, // Object for saving the final results
	const Eigen::MatrixXd& X_in, // Data Matrix
	const Eigen::VectorXi& delay, // Vector indicating the obseravtion index for missing observations for each variable x_n in X
	const int& R, // Number of factors to be estimated by Sparse PCA
	const int& order, // Max order of VAR process in state-equation
	const bool& decorr_errors, // Decorrelate the idiosyncratic errors
	const char* crit, // IC criterion for estimation of the lag order of the VAR process in the state equation
	const double& comp_null, // computational zero
	const bool& check_rank, // check for rank deficiency
	const bool& log, // talk to me
	const double& KFS_conv_crit, //KFS conversion criterion
	const unsigned& fcast_horizon
)
{
	/* Dummies */

	// Integers
	int N = X_in.cols(), T = X_in.rows();

  
  // Matricies
  Eigen::MatrixXd Lambda_hat(N, R), S(N, N), X = X_in, X_PCA = X_in(seq(0, T - delay.maxCoeff() - 1), all);
  Eigen::MatrixXd F_hat(R, X_PCA.rows());
  X.transposeInPlace();
  
  // PCA using the Eigen decomposition
  S = DataHandle::cov(X_PCA);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(S);
  Lambda_hat = eig.eigenvectors().rightCols(R);
  for (int k = 0; k < R; ++k)
  {
    Lambda_hat.col(k).normalize();
  }
  F_hat = (X_PCA * Lambda_hat).transpose();


	// Get order of factor AVR process
	int o = VARorder<BIC>(F_hat, order, comp_null);

	// Calculate inital values for the Kalman filter
	Eigen::MatrixXd F_l = Eigen::MatrixXd::Zero(R * o, T - delay.maxCoeff()), Pt0 = Eigen::MatrixXd::Zero(R * o, R * o), Lambda_l = Eigen::MatrixXd::Zero(N, R * o), Sigma_epsilon = Eigen::MatrixXd::Zero(R * o, R * o),
		Sigma_e = Eigen::MatrixXd::Zero(N, N), Phi = Eigen::MatrixXd::Zero(R * o, R * o);

	for (int oo = 0; oo < o; ++oo)
	{

		F_l.block(R * oo, oo, R, T - delay.maxCoeff() - oo) = F_hat.block(0, 0, R, T - delay.maxCoeff() - oo);

	}

	for (int oo = 1; oo < o; ++oo)
	{

		DataHandle::removeCol(F_l, 0, 0);
		DataHandle::removeCol(X, 0, 0);
		--T;

	}

	{
		// Calculate estimator of the factor VARMA process coefficient matrix Phi
		Eigen::MatrixXd F_t = F_l(Eigen::all, Eigen::seq(1, Eigen::last)).transpose();
		Eigen::MatrixXd F_t_lag = F_l(Eigen::all, Eigen::seq(0, Eigen::last - 1)).transpose();
		Eigen::MatrixXd IdentK = Eigen::MatrixXd::Identity(R * o, R * o);

		Phi = ((F_t_lag.transpose() * F_t_lag).llt().solve(IdentK) * F_t_lag.transpose() * F_t).transpose();
		Phi.bottomLeftCorner(R * (o - 1), R * (o - 1)) = Eigen::MatrixXd::Identity(R * (o - 1), R * (o - 1));
		Phi.bottomRightCorner(R * (o - 1), R) = Eigen::MatrixXd::Zero(R * (o - 1), R);

		Eigen::MatrixXd Phi_T = Phi.transpose();

		// Calculate an estimator for the variance-covariance matrix of the factor model VAR process errors
		for (int t = 0; t < F_t.rows(); ++t)
		{
			Sigma_epsilon.diagonal() += ((F_t.row(t) - F_t_lag.row(t) * Phi_T).transpose() * (F_t.row(t) - F_t_lag.row(t) * Phi_T)).diagonal();
		}
		Sigma_epsilon *= (1. / double(T - 1));

		// Calculate an estimator for the variance-covariance matrix of the idiosyncratic error of the state model
		Lambda_l.topLeftCorner(N, R) = Lambda_hat;

		if (decorr_errors)
		{

			Sigma_e = (X(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1)) - Lambda_l * F_l(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1))) * (X(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1)) - Lambda_l * F_l(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1))).transpose();
		}
		else
		{

			Sigma_e.diagonal() = ((X(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1)) - Lambda_l * F_l(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1))) * (X(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1)) - Lambda_l * F_l(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1))).transpose()).diagonal();

		}
		Sigma_e *= (1. / double(T - delay.maxCoeff()));

		// Calculate initial values for the model uncertainty and factors themselves
		Eigen::MatrixXd F_temp = F_l(Eigen::seq(0, R * o - 1), Eigen::all).transpose();
		//Pt0 = DataHandle::cov(F_temp);
		Pt0.diagonal().setConstant(1000);
	}

	results.F = Eigen::MatrixXd::Zero(R * o, (T + fcast_horizon + 1));
	results.Pt = Eigen::MatrixXd::Zero(R * o, R * o * (T + fcast_horizon + 1));
	results.Wt = Eigen::MatrixXd::Zero(R * o, R * o * T);
	results.Zero_Indeces = Eigen::MatrixXi::Zero(Lambda_hat.rows(), Lambda_hat.cols());

	if (decorr_errors)
	{

		// If the idiosyncratic errors are cross-correlated, transofrm the data
		Eigen::LLT<Eigen::MatrixXd> llt(Sigma_e);

		if (llt.info() != Eigen::Success) {
			Rcpp::Rcerr << "\nWARNING: Decorrelation failed. Using correlated data.\n";
			results.C = Eigen::MatrixXd::Identity(N, N);
		}
		else {
			results.C = llt.matrixL().solve(Eigen::MatrixXd::Identity(N, N));
		}

		if (0 < delay.maxCoeff())
		{

			Eigen::MatrixXd Temp = X;

			for (int t = 0; t < T; ++t)
			{

				Eigen::ArrayXd Ind = (!((Eigen::VectorXi::Constant(N, 1, T) - delay).array() <= t)).cast<double>();

				for (int n = 0; n < N; ++n)
				{
					if (T - delay(n) <= t)
					{
						Temp(n, t) = DBL_MAX;
					}
					else
					{

						Eigen::VectorXd Cn = (results.C.row(n).transpose().array() * Ind).matrix();
						Eigen::VectorXd Xt = X.col(t).array().isNaN().select(0.0, X.col(t)).matrix();

						Temp(n, t) = Cn.dot(Xt);
					}
				}
			}

			X = Temp;

		}
		else
		{
			X = results.C * X;
		}



		Lambda_l.topLeftCorner(N, R) = (results.C * Lambda_l.topLeftCorner(N, R)).eval();


		// Use the diagonalised variance covariance matrix from hereon
		Sigma_e.setZero();
		Sigma_e = (1. / double(T - delay.maxCoeff())) * (X(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1)) - Lambda_l * F_l(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1))) * (X(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1)) - Lambda_l * F_l(Eigen::all, Eigen::seq(0, T - delay.maxCoeff() - 1))).transpose();

	}
	else
	{
		results.C = Eigen::MatrixXd::Identity(N, N);
	}

	/* Filtering */

	Eigen::MatrixXd Vt_inv = Eigen::MatrixXd::Zero(N, T), Kt = Eigen::MatrixXd::Zero(R * o, N * T), e = Eigen::MatrixXd::Zero(N, T), Ft0 = F_l.rowwise().mean();
	UVMVKalmanFilter(results, Vt_inv, Kt, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, Pt0, Ft0, delay, fcast_horizon);

	// Re-try if the filter did obviously not converge
	int tries = 10;
	while (KFS_conv_crit < results.F.cwiseAbs().colwise().mean().mean() && 0 < tries)
	{
		Pt0 *= 1.5;
		UVMVKalmanFilter(results, Vt_inv, Kt, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, Pt0, Ft0, delay, fcast_horizon);
		--tries;

	}
	
	UVMVKalmanSmoother(results, Kt, Vt_inv, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, delay);
	
	results.Lambda_hat = Lambda_l;
	results.order = o;

	return;
}