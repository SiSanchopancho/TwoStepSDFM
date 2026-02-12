# SPDX-License-Identifier: GPL-3.0-or-later
#
#  Copyright \u00A9 2024 Domenic Franjic
#
#  This file is part of TwoStepSDFM.
#
#  TwoStepSDFM is free software: you can redistribute
#  it and/or modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  TwoStepSDFM is distributed in the hope that it
#  will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with TwoStepSDFM. If not, see <https://www.gnu.org/licenses/>.

library(rstudioapi)
library(roxygen2)
library(Rcpp)
library(devtools)
setwd(dirname(getActiveDocumentContext()$path))
rm(list = ls())

# Clean up
unlink("src/*.o")
unlink("src/*.so")
unlink("src/*.dll")
devtools::clean_dll()
devtools::clean_vignettes()

# Compile-test setup (optional)
Rcpp::compileAttributes()
options(pkg.build_extra_flags = FALSE)
devtools::document()
devtools::build_vignettes()
devtools::check(args = "--as-cran")

# Build and install package
devtools::build()
install.packages("../TwoStepSDFM_0.1.5.tar.gz", repos = NULL, type = "source")

