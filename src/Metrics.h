/* -------------------------------------------------------------------------
    This code is part of ElasticScattering.

    Copyright(C) 2022 Stijn Hinlopen

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

#include <windows.h>
#include <vector>

#include "src/sim/es/particle_metrics.h"

struct Metrics {
    ParticleMetrics particle_metrics;
    int real_particles;

    double time_elapsed_lifetimes;
    double time_elapsed_temperatures;

    double avg_particle_lifetime;

    Metrics() {
        real_particles = 0;

        avg_particle_lifetime = 0;

        time_elapsed_lifetimes = 0;
        time_elapsed_temperatures = 0;
    }
};

struct GlobalMetrics
{
    int particles_per_row;
    int phi_steps;
    int unique_impurity_count;
    int cells_per_row;
    double grid_creation_time;
};

struct SampleMetrics
{
    int sample_index;
    int total_particles;
    int total_cells;
    int impurity_count;

    unsigned int seed;
    int total_indexed_impurities;

    bool coherent;

    std::vector<Metrics> iteration_metrics;

    SampleMetrics(int p_sample_index, bool p_coherent, int nlt) {
        sample_index = p_sample_index;
        coherent = p_coherent;

        total_cells = 0;
        total_particles = nlt;
        total_indexed_impurities = 0;

        seed = 0;
        impurity_count = 0;
    }
};
