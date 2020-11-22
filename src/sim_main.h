#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "scattering/ImpurityIndex.h"

#include <string>

void sim_main(const InitParameters& init);
SimulationConfiguration ParseConfig(const std::string file);
std::string GetAvailableDirectory(std::string base);

SigmaResult RunIteration(const std::string output_dir, const int iteration, const ScatteringParameters& sp, const ImpurityIndex& grid);
SimulationResult RunSimulation(SimulationConfiguration& sp);

void LogResult(const SimulationConfiguration& sim_params, const SimulationResult& sr);
void LogImage(const std::string file, const int dim, const double scale, const std::vector<double> data);