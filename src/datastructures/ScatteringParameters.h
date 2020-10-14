#pragma once

typedef struct
{
    int mode;
    int dim; //-
    int impurity_count;
    int integrand_steps;

    int is_clockwise;
    int is_incoherent;
    int is_diag_regions;

    unsigned int impurity_seed;

    double region_size;
    double region_extends;
    double impurity_density;
    double particle_speed;      // v
    double impurity_radius;     // r
    double tau;
    double temperature;

    double alpha;
    double phi;
    double magnetic_field;      // B
    double angular_speed;       // w
} ScatteringParameters;
