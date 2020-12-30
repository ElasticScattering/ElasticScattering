#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "Metrics.h"

#include <string>
#include <fstream>

class Logger {
private:
    static void WriteImageSection(std::vector<unsigned char>& pixels, const std::vector<double>& values, const int dim, const int image_id, const bool colored = false);

    //static void WriteMetricRow(const std::ofstream &file, const std::string metric, float value, int precision);

public:
    static void CreateResultLog(std::string file_path, const SimulationConfiguration& settings, double temperature);
    static void LogResult(const std::string file_path, const DataRow row);
    static void LogImages(const std::string file_path, const int dim, const IterationResult iteration);

    static void CreateMetricsLog(const std::string file_path, const GlobalMetrics& gm);
    static void LogSampleMetrics(const std::string file_path, SampleMetrics metrics);
};