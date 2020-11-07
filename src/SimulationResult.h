#pragma once

#include <vector>

struct SigmaResult {
    double xx;
    double xy;

    SigmaResult() {
        xx = 0;
        xy = 0;
    }
};

struct DataRow {
    double temperature;
    double magnetic_field;

    SigmaResult incoherent;
    SigmaResult coherent;
    double xxd;
};

struct SimulationResult {
    double time_elapsed;

    std::vector<DataRow> results;

    SimulationResult(int size) {
        time_elapsed = 0;
        results.resize(size);
    }
}; 

struct ResultBuffer
{
    std::vector<SigmaResult> intermediate_results;

    ResultBuffer(int size) {
        intermediate_results.resize(size);
    }
};
