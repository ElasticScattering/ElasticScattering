#pragma once

#ifndef DEVICE_PROGRAM
#include <windows.h>
#include <vector>
#endif // !DEVICE_PROGRAM

typedef struct ParticleMetrics {
    int cells_passed;
    int impurities_tested;

    int particles_inside_impurity;
    int particles_escaped;
    int particles_at_bound;

    int max_impurities_tested;
    int max_cells_passed;

#ifndef DEVICE_PROGRAM
    ParticleMetrics() {
        cells_passed = 0;
        impurities_tested = 0;
        particles_inside_impurity = 0;
        particles_escaped = 0;
        particles_at_bound = 0;
        max_impurities_tested = 0;
        max_cells_passed = 0;
    }
#endif
} ParticleMetrics;

#ifndef DEVICE_PROGRAM
typedef struct Metrics {
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
} Metrics;
#endif // !DEVICE_PROGRAM

#ifndef DEVICE_PROGRAM
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
#endif // !DEVICE_PROGRAM