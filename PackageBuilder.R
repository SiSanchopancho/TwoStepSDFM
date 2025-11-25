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

## Set directory auf dieses Skript
setwd(dirname(getActiveDocumentContext()$path))
rm(list = ls())

## Clean up alte Compilate
unlink("src/*.o")
unlink("src/*.so")
unlink("src/*.dll")
devtools::clean_dll()

## Rcpp-Attribute neu erzeugen
Rcpp::compileAttributes()

## Rd + NAMESPACE aus roxygen erzeugen
roxygen2::roxygenise()

## Paket testen (wie CRAN)
devtools::check(args = "--as-cran")

## Build
devtools::build()

# Install
install.packages("../TwoStepSDFM_0.1.4.tar.gz", repos = NULL, type = "source")