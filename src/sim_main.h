#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"

#include <string>

void sim_main(const InitParameters& init);
SimulationConfiguration& ParseConfig(const std::string file);
SimulationResult& RunSimulation(SimulationConfiguration& sp);

void LogResult(const SimulationConfiguration& sim_params, const SimulationResult& sr);
