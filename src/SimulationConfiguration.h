#pragma once

#include "scattering/escl/ScatteringParameters.h"

typedef struct SimulationConfiguration
{
    int number_of_runs;
    int samples_per_run;
    double magnetic_field_min;
    double magnetic_field_max;

    ScatteringParameters scattering_params;
} SimulationConfiguration;
