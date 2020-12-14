#pragma once

#include <windows.h>

class Metrics {
public:
    double time_elapsed;

    int particles_escaped;
    int cells_passed;
    int impurity_intersections;

    int actual_impurity_count;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
};