#pragma once

#include <vector>

struct SimulationResult {
    double time_elapsed;

    std::vector<DataRow> results;

    SimulationResult(int size) {
        time_elapsed = 0;
        results.resize(size);
    }
}; 

struct DataRow {
    double temperature;
    double magnetic_field;
    
    ScatterResult incoherent;
    ScatterResult coherent;
    double xxd;
};

struct ResultBuffer
{
    std::vector<ScatterResult> intermediate_results;

    ResultBuffer(int size) {
        intermediate_results.resize(size);
    }
};

struct ScatterResult {
    double xx;
    double xy;
};
