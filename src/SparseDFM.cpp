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


#include "Internals/SparseDFM.h"
#include <Rcpp.h>

 /* Two-Step Estimator for a LInear Gaussian Sparse Latent (Static) Dynamic Factor Model */
void SparseDFM::SDFMKFS(
	Filtering::KFS_fit& results, // Object for saving the final results
	const Eigen::MatrixXd& X_in, // Data Matrix
	const Eigen::VectorXi& delay, // Vector indicating the obseravtion index for missing observations for each variable x_n in X
	const Eigen::VectorXi& selected, // number of variables ought to be selected for each factor (stopping rule)
	const int& R, // Number of factors to be estimated by Sparse PCA
	const int& order, // Max order of VAR process in state-equation
	const bool& decorr_errors, // Decorrelate the idiosyncratic errors
	const char* crit, // IC criterion for estimation of the lag order of the VAR process in the state equation
	const double& l2, // ridge penalty for the elastic net
	Eigen::VectorXd l1, // lasso penalty of the elastic net (vector of dimensions kx1 since for each respective factor a lasso penalty should be provided) (stopping rule)
	const int& max_iterations, // Maximum number of iterations in the convergence loop in the CD program (currently disabled)
	int steps, // number of steps taken (stopping rule; only for method "LARS")
	const double& comp_null, // computational zero
	const bool& check_rank, // check for rank deficiency
	const double& conv_crit, // conversion criterion for the sparse pca algorithm
	const double& conv_threshold, // conversion criterion for the conversion loop in the CD program (currently disabled)
	const bool& log, // talk to me
	const double& KFS_conv_crit //KFS conversion criterion
)
{
	/* Dummies */

	// Integers
	int N = X_in.cols(), T = X_in.rows();

	// Matricies
	Eigen::MatrixXd X = X_in;

	X.transposeInPlace();

	// Sparse PCA
	SparsePCA::SPC_fit fit;
	{
		SparsePCA::SparsePC(fit, X_in(Eigen::seq(0, T - delay.maxCoeff() - 1), Eigen::all), selected, R, l2, l1, max_iterations, steps, comp_null,
			check_rank, conv_crit, conv_threshold);
	}


	// Get order of factor AVR process
	int o = VARorder<BIC>(fit.F_hat, order, comp_null);

	// Calculate inital values for the Kalman filter
	Eigen::MatrixXd F_l = Eigen::MatrixXd::Zero(R * o, T - delay.maxCoeff()), Pt0 = Eigen::MatrixXd::Zero(R * o, R * o), Lambda_l = Eigen::MatrixXd::Zero(N, R * o), Sigma_epsilon = Eigen::MatrixXd::Zero(R * o, R * o),
		Sigma_e = Eigen::MatrixXd::Zero(N, N), Phi = Eigen::MatrixXd::Zero(R * o, R * o);

	for (int oo = 0; oo < o; ++oo)
	{

		F_l.block(R * oo, oo, R, T - delay.maxCoeff() - oo) = fit.F_hat.block(0, 0, R, T - delay.maxCoeff() - oo);

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
		Lambda_l.topLeftCorner(N, R) = fit.Lambda_hat;

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
		Pt0 = DataHandle::cov(F_temp);
	}

	results.F = Eigen::MatrixXd::Zero(R * o, (T + 0 + 1));
	results.Pt = Eigen::MatrixXd::Zero(R * o, R * o * (T + 0 + 1));
	results.Wt = Eigen::MatrixXd::Zero(R * o, R * o * T);
	results.Zero_Indeces = Eigen::MatrixXi::Zero(fit.Lambda_hat.rows(), fit.Lambda_hat.cols());

	if (decorr_errors)
	{

		// Save the zero indices for later
		results.Zero_Indeces = fit.Lambda_hat.unaryExpr([comp_null](double a) { return (a == 0) ? 1 : 0; });


		// If the idiosyncratic errors are cross-correlated, transofrm the data
		Eigen::LLT<Eigen::MatrixXd> llt(Sigma_e);

		results.C = llt.matrixL().solve(Eigen::MatrixXd::Identity(N, N));

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
						Temp(n, t) = NAN;
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



		Lambda_l.topLeftCorner(N, R) = results.C * Lambda_l.topLeftCorner(N, R);


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

	UVMVKalmanFilter(results, Vt_inv, Kt, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, Pt0, Ft0, delay, 0);
	UVMVKalmanSmoother(results, Kt, Vt_inv, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, delay);

	// Re-try if the filter did obviously not converge
	int tries = 10;
	while (KFS_conv_crit < results.F.cwiseAbs().colwise().mean().mean() && 0 < tries)
	{
		Pt0 *= 0.5;
		UVMVKalmanFilter(results, Vt_inv, Kt, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, Pt0, Ft0, delay, 0);
		UVMVKalmanSmoother(results, Kt, Vt_inv, e, N, T, R * o, X, Lambda_l, Phi, Sigma_e, Sigma_epsilon, delay);
		--tries;

	}

	results.Lambda_hat = Lambda_l;
	results.order = o;

	return;
}