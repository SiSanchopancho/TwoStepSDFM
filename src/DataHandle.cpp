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


#include "Internals/DataHandle.h"

/* Compute the variance co-variance matrix */
Eigen::MatrixXd DataHandle::cov(const Eigen::MatrixXd& X_in)
{
    Eigen::MatrixXd X = X_in;
    Eigen::VectorXd mu = X.colwise().mean();
    for (int n = 0; n < X.cols(); ++n)
    {
        for (int t = 0; t < X.rows(); ++t)
        {
            X(t, n) -= mu(n, 0);
        }
    }
    return (X.transpose() * X) / double(X.rows() - 1);
};

/* Remove row from matrix in-place */
void DataHandle::removeRow(Eigen::MatrixXd& X, const int& t, const bool& conservative)
{
    int T = static_cast<int>(X.rows()) - 1;
    int N = static_cast<int>(X.cols());
    if (t < T)
    {
        X.block(t, 0, T - t, N) = X.block(t + 1, 0, T - t, N).eval();
    }
    if (conservative)
    {
        X.row(T).setZero();
    }
    else
    {
        X.conservativeResize(T, N);
    }

};

/* Remove row from matrix in-place*/
void DataHandle::removeCol(Eigen::MatrixXd& X, const int& n, const bool& conservative)
{
    int T = static_cast<int>(X.rows());
    int N = static_cast<int>(X.cols()) - 1;
    if (n < N)
    {
        X.block(0, n, T, N - n) = X.block(0, n + 1, T, N - n).eval();
    }
    if (conservative)
    {
        X.col(N).setZero();
    }
    else
    {
        X.conservativeResize(T, N);
    }
};
