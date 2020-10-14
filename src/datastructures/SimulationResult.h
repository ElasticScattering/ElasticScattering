#pragma once

#include <vector>

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
        time_elapsed = 0;
        results.resize(size);
    }
};