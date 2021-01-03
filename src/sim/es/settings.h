#pragma once
#ifndef DEVICE_PROGRAM
    #include "v2.h"
#endif

typedef struct ImpuritySettings
{
    double spawn_region_start;
    double spawn_region_size;

    double cell_size;
    int    cells_per_row;

    double impurity_radius;
} ImpuritySettings;

typedef struct Range
{
    double start;
    double step_size;
    double count;
} Range;

typedef struct ParticleSettings
{
    double particle_speed;
    double angular_speed;
    double alpha;

    int is_clockwise;
    int is_coherent;

    double phi_step_size;
    double phi_start;
} ParticleSettings;

typedef struct SimulationSettings
{
    int particles_per_quadrant;
    int particles_per_position;
    int particles_per_row;
    int total_particles;

    int positions_per_row;
    int total_positions;

    double signed_angular_speed; //?

    double integrand_angle_area;
    double phi_integrand_factor;

    double region_size;
    double region_extended_area;
    
    double distance_between_positions;
    double2 small_offset;

    double coherent_tau;
} SimulationSettings;
