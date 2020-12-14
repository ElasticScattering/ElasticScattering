#pragma once

#include <vector>

struct SigmaResult {
    std::vector<double> xx_buffer, xy_buffer;
};

struct Sigma {
    double xx, xy;

    Sigma() { xx = 0; xy = 0; }

    Sigma(double _xx, double _xy) { xx = _xx; xy = _xy; }

    Sigma operator+(const Sigma& a) const
    {
        return Sigma(a.xx + xx, a.xy + xy);
    }

    Sigma& operator=(const Sigma& a)
    {
        xx = a.xx;
        xy = a.xy;
        return *this;
    }

    Sigma& operator+=(const Sigma& a)
    {
        xx += a.xx;
        xy += a.xy;
        return *this;
    }
};

struct DataRow {
    double temperature;
    double magnetic_field;
    
    Sigma incoherent;
    Sigma coherent;
    double xxd;

    DataRow() {}

    DataRow(Sigma _coherent, Sigma _incoherent, int samples)
    {
        double s = (double)(samples);
        coherent.xx = _coherent.xx / s;
        coherent.xy = _coherent.xy / s;

        incoherent.xx = _incoherent.xx / s;
        incoherent.xy = _incoherent.xy / s;
    }
};

struct IterationResult {
    std::vector<double> particle_lifetimes;
    SigmaResult sigmas;
    Sigma result;

    IterationResult() {}
};


class SampleResult {
public:
    void AddSample();
};

class SimulationResult {
    double time_elapsed;
    //SimulationConfiguration cfg;
    std::vector<SampleResult> samples;
};