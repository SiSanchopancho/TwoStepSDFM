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

#ifndef SDFM
#define SDFM

#define _USE_MATH_DEFINES

// Including external libraries
#include <RcppEigen.h>
#include <stdlib.h>
#include <math.h>

// Including internal libraries
#include "Filtering.h"
#include "SparsePCA.h"
#include "Orders.h"

namespace SparseDFM {

    extern
        void SDFMKFS(
            Filtering::KFS_fit& results,
            const Eigen::MatrixXd& X_in,
            const Eigen::VectorXi& delay,
            const Eigen::VectorXi& selected,
            const int& R = 1,
            const int& order = 10,
            const bool& decorr_errors = 0,
            const char* crit = "BIC",
            const double& l2 = 0.4,
            Eigen::VectorXd l1 = INT_MIN * Eigen::VectorXd::Zero(1),
            const int& max_iterations = 1000,
            int steps = INT_MIN,
            const double& comp_null = 10e-15,
            const bool& check_rank = 0,
            const double& conv_crit = 0.0001,
            const double& conv_threshold = 10e-7,
            const bool& log = 0,
            const double& KFS_conv_crit = 10e+2
        );

};
#endif /* defined(SDFM) */

