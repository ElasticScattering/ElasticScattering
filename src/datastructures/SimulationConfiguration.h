#pragma once

#include "ScatteringParameters.h"

typedef struct SimulationConfiguration
{
    ScatteringParameters scattering_params;
    int runs;
    int samples_per_run;
    double magnetic_field_min;
    double magnetic_field_max;

} SimulationConfiguration;
