#pragma once
#ifndef DEVICE_PROGRAM
    #include "v2.h"
#endif

typedef struct ImpuritySettings
{
    double2 impurity_spawn_range;

    double impurity_radius;

    double cell_size;
    int cells_per_row;
    int impurity_count;
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
