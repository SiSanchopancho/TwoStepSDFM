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
#include <ostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <chrono>
#include <fstream>
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h> // For _mkdir on Windows
#define mkdir _mkdir // Define mkdir to use _mkdir on Windows
#else
#include <unistd.h> // For access and mkdir on POSIX systems
#endif

using namespace Eigen;

namespace DataHandle {

    Eigen::MatrixXd cov(const Eigen::MatrixXd& X_in);
    void removeRow(Eigen::MatrixXd& X, const int& t, const bool& conservative = true);
    void removeCol(Eigen::MatrixXd& X, const int& n, const bool& conservative = true);

};

#ifdef _MSC_VER // Check if the compiler is MSVC
#pragma warning(disable : 4996) // Disable MSVC warning due to strcpy
#endif

/** Save matrix to CSV */
template <typename Derived>
void dataSaveCSV(const Eigen::MatrixBase<Derived>& matrix, const std::string& name, const std::string& file_format = ".txt")
{

	std::string out_name = name + file_format; // Create the complete file name
	const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n"); // Open a csv IO-format
	std::ofstream file;
	file.open(out_name, std::ofstream::out | std::ofstream::trunc);
	if (file.is_open())
	{
		file << matrix.format(CSVFormat);
	}
	file.close();
};

/** Load matrix from CSV */
inline Eigen::MatrixXd dataLoadCSV(const std::string& name)
{
	std::vector<double> entry;
	std::string rows, element;
	std::ifstream matrixDataFile(name);
	unsigned t = 0;
	while (getline(matrixDataFile, rows))
	{
		std::stringstream matrixRowStringStream(rows);
		while (getline(matrixRowStringStream, element, ','))
		{
			entry.push_back(stod(element));
		}
		++t;
	}

	// Map std::vector data pointer to matrix and return
	return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(entry.data(), t, entry.size() / t);
};

/** Load matrix from CSV that contains a single column of data */
inline Eigen::MatrixXd dataLoadCSV2(const std::string& name)
{
	std::vector<double> entry;
	std::string rows;
	std::ifstream matrixDataFile(name);
	while (getline(matrixDataFile, rows))
	{
		entry.push_back(stod(rows));
	}
	return Map<Matrix<double, Dynamic, Dynamic, ColMajor>>(entry.data(), entry.size(), 1);
};

/** Load matrix from CSV that conatins date data */
inline Eigen::MatrixXd dataLoadCSVDate(const std::string& name)
{

	std::vector<std::string> dateStrings;
	std::string rows;
	std::ifstream matrixDataFile(name);

	// Read the data by element since the .csv consists of a single line of data
	while (getline(matrixDataFile, rows))
	{
		dateStrings.push_back(rows);
	}

	std::vector<std::tm> dates;
	for (const auto& dateString : dateStrings) {

		std::tm date = {};
		std::istringstream ss(dateString);
		ss >> std::get_time(&date, "%Y-%m-%d");
		if (ss.fail()) {
			throw std::invalid_argument("Error parsing date.");
		}
		dates.push_back(date);

	}

	// Convert std::vector<std::tm> to Eigen::MatrixXd
	Eigen::MatrixXd result(dates.size(), 3);
	for (size_t i = 0; i < dates.size(); ++i)
	{
		result(i, 0) = dates[i].tm_year + 1900;  // Years since 1900
		result(i, 1) = dates[i].tm_mon + 1;      // Months are 0-based
		result(i, 2) = dates[i].tm_mday;    // Months are 0-based
	}

	return result;

}

/** Check whether a given directory exists */
inline bool folderExists(const std::string& folderPath) {
#ifdef _WIN32 // Existence check on Windows machiens
	struct _stat info;
	return _stat(folderPath.c_str(), &info) == 0 && (info.st_mode & _S_IFDIR);
#else // Existence check on Linux/Ubuntu machiens
	struct stat info;
	return stat(folderPath.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
#endif
}

/** Create a folder at the given path */
inline bool createFolder(const std::string& folderPath) {
#ifdef _WIN32 // Create folder on Windows machines
	return _mkdir(folderPath.c_str()) == 0;
#else // Create folder on Linux/Ubutnu machines
	return mkdir(folderPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0;
#endif
}

/** Ensure that a folder exists, create it if it does not */
inline bool ensureFolderExists(const std::string& folderName) {

	if (!folderExists(folderName)) {
		if (createFolder(folderName)) {
			Rcpp::Rcout << "\nFolder '" << folderName << "' created successfully." << std::endl;
			return true;
		}
		else {
			throw std::invalid_argument("Folder cannot be created!");
		}
	}
	else {
		Rcpp::Rcout << "\nFolder '" << folderName << "' already exists." << std::endl;
		return true;
	}
}

/** Get a string of the current date formatted as dd-mm-yyyy */
inline std::string getFormattedDate() {

	auto now = std::chrono::system_clock::now();
	std::time_t t = std::chrono::system_clock::to_time_t(now);
	std::tm tm_struct;

#ifdef _WIN32 // Retrieve the current time on Windows
	localtime_s(&tm_struct, &t);
#else // Retrieve the current time on Ubuntu/Linux
	localtime_r(&t, &tm_struct);
#endif

	std::ostringstream oss;
	oss << std::put_time(&tm_struct, "%d-%m-%Y");  // Format as dd-mm-yyyy

	return oss.str();
}

/** Create a sequence of logarithmically spaced values */
inline Eigen::VectorXd logSeq(const double& start, const double& end, const unsigned& length)
{

	Eigen::VectorXd lseq = Eigen::VectorXd::Zero(length);
	double curr = start;
	double step = std::pow((end / start), double(1.0 / double(length - 1)));

	for (int i = 0; i < length; ++i)
	{
		lseq(i) = curr;
		curr *= step;
	}

	return lseq;
}


#endif /* defined(DATA_HANDLING) */
