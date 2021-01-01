#pragma once

#ifndef DEVICE_PROGRAM
#include <windows.h>
#include <vector>
#endif // !DEVICE_PROGRAM

struct Metrics {
    int cells_passed;
    int impurities_tested;
    int particles_inside_impurity;
    int particles_escaped;

    int max_impurities_tested;
    int max_cells_passed;

#ifndef DEVICE_PROGRAM
    int real_lifetimes;

    double time_elapsed_lifetimes;
    double time_elapsed_temperatures;

    double avg_particle_lifetime;


    Metrics() {
        cells_passed = 0;
        impurities_tested = 0;
        particles_inside_impurity = 0;
        real_lifetimes = 0;
        particles_escaped = 0;

        avg_particle_lifetime = 0;

        time_elapsed_lifetimes = 0;
        time_elapsed_temperatures = 0;

        max_cells_passed = 0;
        max_impurities_tested = 0;
    }
#endif // !DEVICE_PROGRAM
};

#ifndef DEVICE_PROGRAM
struct GlobalMetrics
{
    int particles_per_row;
    int phi_values;
    int unique_impurity_count;
    int cells_per_row;
    double grid_creation_time;
};

struct SampleMetrics
{
    int sample_index;
    int nlifetimes;
    int cells_per_row;
    int impurity_count;

    int seed;
    int total_indexed_impurities;

    bool coherent;

    std::vector<Metrics> iteration_metrics;

    SampleMetrics(int p_sample_index, bool p_coherent, int N, int nlt) {
        sample_index = p_sample_index;
        coherent = p_coherent;
        iteration_metrics.resize(N);

        nlifetimes = nlt;
    }
};
#endif // !DEVICE_PROGRAM