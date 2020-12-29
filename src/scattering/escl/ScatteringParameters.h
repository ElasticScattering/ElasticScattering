#pragma once
#ifndef DEVICE_PROGRAM
    #include "v2.h"
#endif

typedef struct ImpuritySettings
{
    double spawn_region_start;
    double spawn_region_size;

    double cell_size;
    int cells_per_row;

    double impurity_radius;
} ImpuritySettings;

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
