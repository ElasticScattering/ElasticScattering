#pragma once

#include <vector>
#include "src/SimulationParameters.h"

struct SimulationResult {
    double time_elapsed;

    bool x_is_temperature;

    std::vector<double> xs;
    std::vector<double> xs_temperature;

    std::vector<double> results_xx;
    std::vector<double> results_xy;
    std::vector<double> results_xxi;
    std::vector<double> results_xyi;

    std::vector<double> delta_xxi;

    SimulationResult(int size) {
        xs.resize(size);
        xs_temperature.resize(size);
        results_xx.resize(size);
        results_xxi.resize(size);
        delta_xxi.resize(size);
        results_xy.resize(size);
        results_xyi.resize(size);
        x_is_temperature = false;
    }
};