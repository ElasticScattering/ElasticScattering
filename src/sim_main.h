#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "scattering/ElasticScattering.h"

int sim_main(const InitParameters& init);
SimulationResult& RunSimulation(ElasticScattering& es, SimulationConfiguration& sp);
void PrintInfo(const SimulationConfiguration& sp, int count);
void LogResult(const SimulationConfiguration& sim_params, const SimulationResult& sr);
