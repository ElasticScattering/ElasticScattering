/* -------------------------------------------------------------------------
	This code is part of ElasticScattering.

	Copyright(C) 2022 Stijn Hinlopen

	This program is free software : you can redistribute it and /or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "Metrics.h"

#include <string>
#include <fstream>

class Logger {
private:
    static const int result_precision = 6;
    static const int rw = result_precision + 10;
    static void WriteImageSection(std::vector<unsigned char>& pixels, const std::vector<double>& values, const int dim, const int image_id, const bool colored = false);

public:
    static void CreateResultLog(const std::string file_path, const SimulationConfiguration& settings, double temperature);
    static void LogResult(const std::string file_path, const DataRow& row);

    static void CreateSampleResultLog(const std::string file_path);
    static void LogSampleResults(const std::string file_path, const SampleResult& coherent, const SampleResult& incoherent);

    static void CreateSampleMetricsLog(const std::string file_path, const GlobalMetrics& gm);
    static void LogSampleMetrics(const std::string file_path, const SampleMetrics& metrics);

    static void LogImages(const std::string file_path, const int dim, const IterationResult& iteration);
};