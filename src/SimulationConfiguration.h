#pragma once

#include "Settings.h"

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
    std::vector<double> magnetic_fields;
    std::vector<double> temperatures;

    int particles_per_row;
    int quadrant_integral_steps;

    Settings settings;
} SimulationConfiguration;
