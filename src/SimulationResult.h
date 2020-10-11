#pragma once

#include <vector>
#include "src/SimulationParameters.h"


struct Result {
    double x;
    double xx;
    double xxi;
    double xxd;
    double xy;
    double xyi;
};

struct SimulationResult {
    double time_elapsed;

    std::vector<Result> results;

    SimulationResult(int size) {
        results.resize(size);
    }
};