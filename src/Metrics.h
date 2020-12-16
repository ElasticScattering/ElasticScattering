#pragma once

#include <windows.h>

struct Metrics {
    int cells_passed;
    int impurity_intersections;
    
    int particles_inside_impurity;
    double pct_particles_inside_impurity;

    int particles_escaped;
    double pct_particles_escaped;

    double prt_cells_passed;
    double prt_impurity_intersections;

    double pct_prt_cells_passed;
    double pct_prt_impurity_intersections;

    double time_elapsed_lifetimes;
    double time_elapsed_temperatures;
};

class SampleMetrics {
public:
    Metrics coherent, incoherent;

    int actual_impurity_count;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
};