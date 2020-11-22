#pragma once

#include "scattering/escl/ScatteringParameters.h"

enum ProgramMode {
    Test,
    Simulation
};

typedef struct
{
    ProgramMode mode;
    bool use_gpu;
    bool dont_show_info;
    bool write_images;
} InitParameters;

typedef struct SimulationConfiguration
{
    int number_of_runs;
    int samples_per_run;
    double magnetic_field_min;
    double magnetic_field_max;

    ScatteringParameters scattering_params;
} SimulationConfiguration;
