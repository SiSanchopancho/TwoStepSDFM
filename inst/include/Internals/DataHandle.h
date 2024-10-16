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

#ifndef DATA_HANDLING
#define DATA_HANDLING

#ifdef _MSC_VER // Check if the compiler is MSVC
#pragma warning(disable : 4996) // Disable MSVC warning due to strcpy
#endif

// Including external libraries
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h> // For _mkdir on Windows
#define mkdir _mkdir // Define mkdir to use _mkdir on Windows
#else
#include <unistd.h> // For access and mkdir on POSIX systems
#endif

namespace DataHandle {

    Eigen::MatrixXd cov(const Eigen::MatrixXd& X_in);
    void removeRow(Eigen::MatrixXd& X, const int& t, const bool& conservative = true);
    void removeCol(Eigen::MatrixXd& X, const int& n, const bool& conservative = true);

};
#endif /* defined(DATA_HANDLING) */
