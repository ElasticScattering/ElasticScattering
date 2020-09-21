#pragma once

#include <vector>
#include "src/SimulationParameters.h"

struct SimulationResult {
    int n_runs;
    int iterations_per_run;
    double time_elapsed;

    bool x_is_temperature;

    std::vector<double> xs;
    std::vector<double> xs_temperature;

    std::vector<double> results_xx;
    std::vector<double> results_xy;
    std::vector<double> results_xxi;
    std::vector<double> results_xyi;

    std::vector<double> delta_xxi;
};