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
    double min, max, step_size;
    int n;
} Range;

typedef struct SimulationConfiguration
{
    int num_samples;
    Range magnetic_field_range;
    std::vector<double> temperatures;

    std::string base_output_directory;
    std::string current_intermediates_dir;
    std::string current_metrics_directory;

    ScatteringParameters scattering_params;
} SimulationConfiguration;
