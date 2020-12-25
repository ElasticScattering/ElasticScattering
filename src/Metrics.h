#pragma once

#ifndef DEVICE_PROGRAM
#include <windows.h>
#endif // !DEVICE_PROGRAM

#ifndef DEVICE_PROGRAM
struct GlobalMetrics
{
    bool coherent;

    int particles_per_row;
    int phi_values;
    int unique_impurity_count;
    int additional_impurities;
    int cells_per_row;

    double grid_time_elapsed;
    double avg_impurities_in_cell;
    double avg_impuritiies_in_cell_overlapping;
};
#endif // !DEVICE_PROGRAM

struct Metrics {
    int cells_passed;
    int impurity_intersections;
    int particles_inside_impurity;
    int particles_escaped;

#ifndef DEVICE_PROGRAM
    int magnetic_field_index;

    double mln_cells_passed;
    double mln_impurity_intersections;
    double pct_particles_inside_impurity;

    double prt_cells_passed;
    double prt_impurity_intersections;

    double pct_prt_cells_passed;
    double pct_prt_impurity_intersections;

    double time_elapsed_lifetimes;
    double time_elapsed_temperatures;

    Metrics(int mf) {
        magnetic_field_index = mf;

        cells_passed = 0;
        impurity_intersections = 0;
        particles_inside_impurity = 0;
        particles_escaped = 0;

        pct_particles_inside_impurity = 0;
        prt_cells_passed = 0;
        prt_impurity_intersections = 0;
        pct_prt_cells_passed = 0;
        pct_prt_impurity_intersections = 0;

        time_elapsed_lifetimes = 0;
        time_elapsed_temperatures = 0;
    }
#endif // !DEVICE_PROGRAM
};
