#pragma once

#ifndef DEVICE_PROGRAM
#include <windows.h>
#endif // !DEVICE_PROGRAM

#ifndef DEVICE_PROGRAM
struct GlobalMetrics
{
    int particles_per_row;
    int phi_values;
    int unique_impurity_count;
    int additional_impurities;
    int cells_per_row;

    double grid_time_elapsed;
    double avg_impurities_in_cell;
    double avg_impurities_in_cell_overlapping;
};
#endif // !DEVICE_PROGRAM

struct Metrics {
    int cells_passed;
    int impurity_intersections;
    int particles_inside_impurity;
    int particles_escaped;

#ifndef DEVICE_PROGRAM
    int magnetic_field_index;
    int nlifetimes;
    int cells_per_row;
    int impurity_count;

    double time_elapsed_lifetimes;
    double time_elapsed_temperatures;

    double avg_particle_lifetime;

    Metrics(int mf, int nlt, int cells, int imp_count) {
        magnetic_field_index      = mf;
        nlifetimes                = nlt;
        cells_per_row             = cells;
        impurity_count            = imp_count;

        cells_passed              = 0;
        impurity_intersections    = 0;
        particles_inside_impurity = 0;
        particles_escaped         = 0;

        avg_particle_lifetime     = 0;

        time_elapsed_lifetimes    = 0;
        time_elapsed_temperatures = 0;
    }
#endif // !DEVICE_PROGRAM
};
