#pragma once

typedef struct Settings
{
    double tau;
    double alpha;
    bool is_clockwise;
    double particle_speed;

    double region_size;
    double region_extends;
    double impurity_density;
    double impurity_radius;
    int target_cell_population;
} Settings;
