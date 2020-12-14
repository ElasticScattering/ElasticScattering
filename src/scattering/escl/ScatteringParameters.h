#pragma once
#ifndef DEVICE_PROGRAM
    #include "v2.h"
#endif

typedef struct ScatteringParameters
{
    int dim;
    int integrand_steps; 
    int values_per_particle;
    int impurity_count;

    double magnetic_field;
    double tau;
    double temperature;
    double default_max_lifetime;

    double particle_speed;
    double angular_speed;
    double alpha;

    int is_clockwise;
    int is_incoherent;
    
    double integrand_step_size;
    double integrand_start_angle;
    double integrand_angle_area;

    double region_size;
    double region_extends;
    double impurity_density;
    double impurity_radius;

    int max_expected_impurities_in_cell;
    int cells_per_row;
    double cell_size;

    double2 impurity_spawn_range;
} ScatteringParameters;
