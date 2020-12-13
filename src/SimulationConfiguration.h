#pragma once

#include "scattering/escl/ScatteringParameters.h"

#include <vector>
#include <string>

enum class ProgramMode {
    Test,
    Simulation
};

typedef struct
{
    ProgramMode mode;
    bool use_gpu;
    bool dont_show_info;
    bool write_images;
    std::string config_file;
} InitParameters;

typedef struct Range
{
    double min, max;
    int n;
    int step_size;
} Range;

typedef struct SimulationConfiguration
{
    int samples_per_run;
    Range magnetic_field_range;
    std::vector<double> temperatures;

    std::string output_directory;
    std::string intermediates_directory;

    ScatteringParameters scattering_params;
} SimulationConfiguration;
