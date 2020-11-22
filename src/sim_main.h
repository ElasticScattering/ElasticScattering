#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "scattering/ElasticScattering.h"

#include <string>

int sim_main(const InitParameters& init);
SimulationResult& RunSimulation(ElasticScattering& es, SimulationConfiguration& sp);
SimulationConfiguration& ParseConfig(const std::string file);

void LogResult(const SimulationConfiguration& sim_params, const SimulationResult& sr);
