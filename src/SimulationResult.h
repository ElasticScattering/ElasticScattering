/* -------------------------------------------------------------------------
    This code is part of ElasticScattering.

    Copyright(C) 2022 Elastic Scattering developers

    This program is free software : you can redistribute it and /or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#pragma once

#include <vector>

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

struct SampleResult {
    std::vector<std::vector<Sigma>> results;

    SampleResult() {}

    SampleResult(int T, int N)
    {
        results.resize(N);
        for (int i = 0; i < N; i++)
            results[i].resize(T);
    }
};

struct SimulationResult {
    std::vector<SampleResult> coherent;
    std::vector<SampleResult> incoherent;

    SimulationResult(const int S)
    {
        coherent.reserve(S);
        incoherent.reserve(S);
    }
};

struct SigmaResult {
    std::vector<double> xx_buffer, xy_buffer;
};

typedef struct IterationResult {
    std::vector<double> particle_lifetimes;
    SigmaResult sigmas;
    Sigma result;

    IterationResult() {}
} IterationResult;


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

    DataRow(double temp, double mf, Sigma _coherent, Sigma _incoherent)
    {
        temperature = temp;
        magnetic_field = mf;
        coherent = _coherent;
        incoherent = _incoherent;
    }
};
