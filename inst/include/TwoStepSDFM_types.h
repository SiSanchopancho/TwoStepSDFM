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

#ifndef KFSMLE_TYPES
#define KFSMLE_TYPES

// Include necessary Rcpp and Eigen libraries
#include <RcppCommon.h>
#include <Rcpp.h>
#include <Eigen/Eigen>

// Results structure to store the outcomes of the optimization process
struct Results {
    double min_val;           // Minimum value found by the optimization
    Eigen::VectorXd estimate; // Estimated parameters
    Eigen::MatrixXd hessian;  // Hessian matrix at the optimum
};

// Rcpp namespace to enable seamless integration with R
namespace Rcpp {
    // Specialize the wrap function for the Results struct
    template <>
    SEXP wrap(const Results& x);
}

#endif /* defined(KFSMLE_TYPES) */
