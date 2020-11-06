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
    
    SigmaResult incoherent;
    SigmaResult coherent;
    double xxd;
};

struct ResultBuffer
{
    std::vector<SigmaResult> intermediate_results;

    ResultBuffer(int size) {
        intermediate_results.resize(size);
    }
};

struct SigmaResult {
    double xx;
    double xy;
};
