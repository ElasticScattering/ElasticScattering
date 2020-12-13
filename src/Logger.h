#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"

#include <string>

class Logger {
private:
    static void WriteImageSection(std::vector<unsigned char>& pixels, const std::vector<double>& values, const int dim, const int image_id, const bool colored = false);

public:
    static void CreateLog(const std::string result_file, const SimulationConfiguration& cfg, double temperature);
    static void FinishLog(const std::string file_path, const double time_elapsed);

    static void LogResult(const std::string file_path, const DataRow row);
    static void LogImages(const std::string file, const int dim, const double scale, const IterationResult iteration);
};