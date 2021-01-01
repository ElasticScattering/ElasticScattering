#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "Metrics.h"

#include <string>
#include <fstream>

class Logger {
private:
    static void WriteImageSection(std::vector<unsigned char>& pixels, const std::vector<double>& values, const int dim, const int image_id, const bool colored = false);

public:
    static void CreateResultLog(const std::string file_path, const SimulationConfiguration& settings, double temperature);
    static void LogResult(const std::string file_path, const DataRow& row);

    static void CreateSampleResultLog(const std::string file_path, const SimulationConfiguration& settings);
    static void LogSampleResults(const std::string file_path, const SampleResult coherent, const SampleResult incoherent);

    static void CreateMetricsLog(const std::string file_path, const GlobalMetrics& gm);
    static void LogSampleMetrics(const std::string file_path, const SampleMetrics& metrics);

    static void LogImages(const std::string file_path, const int dim, const IterationResult& iteration);
};