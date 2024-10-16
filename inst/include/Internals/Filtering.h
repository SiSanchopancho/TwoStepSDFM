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


#pragma once

#ifndef FILTERING
#define FILTERING

// Including external libraries
#include <Eigen/Eigen>
#include <stdlib.h>
#include <math.h>

namespace Filtering {

    class KFS_fit {
    public:
        Eigen::MatrixXd Lambda_hat; // Matrix of estimated factor loadings
        Eigen::MatrixXi Zero_Indeces; // Matrix of estimated factor loadings
        Eigen::MatrixXd Pt; // Conditional variance of factors at t+1 given observations up to t
        Eigen::MatrixXd F; // Smoothed estimates of latent factors
        Eigen::MatrixXd Wt; // Conditional variance of factors at time t given all observations
        Eigen::MatrixXd C; // Inverse of the Lower Cholesky-Deco of Var-Cov
        int order; // Order of VAR process for factors in the transition equation
        bool conv; // Whether the filter converged

        // Default Constructor
        KFS_fit() : order(0), conv(true) {
            Lambda_hat = Eigen::MatrixXd::Zero(0, 0);
            Zero_Indeces = Eigen::MatrixXi::Zero(0, 0);
            Pt = Eigen::MatrixXd::Zero(0, 0);
            F = Eigen::MatrixXd::Zero(0, 0);
            Wt = Eigen::MatrixXd::Zero(0, 0);
            C = Eigen::MatrixXd::Zero(0, 0);
        }

        // Constructor
        KFS_fit(int R, int T, int N) : order(0), conv(false) {
            Lambda_hat = Eigen::MatrixXd::Zero(N, R);
            Zero_Indeces = Eigen::MatrixXi::Zero(N, R);
            Pt = Eigen::MatrixXd::Zero(R, R);
            F = Eigen::MatrixXd::Zero(R, T);
            Wt = Eigen::MatrixXd::Zero(R, R);
            C = Eigen::MatrixXd::Zero(N, N);
        }

        // Copy constructor
        KFS_fit(const KFS_fit& other) = default;

        // Move constructor
        KFS_fit(KFS_fit&& other) noexcept = default;

        // Copy assignment operator
        KFS_fit& operator=(const KFS_fit& other) = default;

        // Move assignment operator
        KFS_fit& operator=(KFS_fit&& other) noexcept = default;

        // Destructor
        ~KFS_fit() = default;
    };


    extern

        void UVMVKalmanFilter(
            Filtering::KFS_fit& results,
            Eigen::MatrixXd& Vt_inv,
            Eigen::MatrixXd& Kt,
            Eigen::MatrixXd& e,
            const int& N,
            const int& T,
            const int& K,
            const Eigen::MatrixXd& X_in,
            const Eigen::MatrixXd& Lambda,
            const Eigen::MatrixXd& Phi,
            const Eigen::MatrixXd& Sigma_e,
            const Eigen::MatrixXd& Sigma_epsilon,
            const Eigen::MatrixXd& Pt0,
            const Eigen::VectorXd& Ft0,
            const Eigen::VectorXi& missings,
            const int& h = 0
        );

    void UVMVKalmanSmoother(
        Filtering::KFS_fit& results,
        Eigen::MatrixXd& Kt,
        Eigen::MatrixXd& Vt_inv,
        Eigen::MatrixXd& e,
        const int& N,
        const int& T,
        const int& K,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Lambda,
        const Eigen::MatrixXd& Phi,
        const Eigen::MatrixXd& Sigma_e,
        const Eigen::MatrixXd& Sigma_epsilon,
        const Eigen::VectorXi& missings
    );
};
#endif /* defined(FILTERING) */

