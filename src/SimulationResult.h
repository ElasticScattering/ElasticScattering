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

    DataRow() {}

    DataRow(SigmaResult _coherent, SigmaResult _incoherent, int samples)
    {
        double s = (double)(samples);
        coherent.xx = _coherent.xx / s;
        coherent.xy = _coherent.xy / s;

        incoherent.xx = _incoherent.xx / s;
        incoherent.xy = _incoherent.xy / s;
    }
};

struct SigmaBuffer {
    std::vector<double> sigma_xx, sigma_xy;

    SigmaBuffer() {}
    SigmaBuffer(std::vector<double>& xx, std::vector<double>& xy) {
        sigma_xx = xx;
        sigma_xy = xy;
    }
};

struct IterationResult {
    SigmaBuffer sigma;
    std::vector<double> lifetimes;
    SigmaResult result;

    IterationResult() {}
};
