#pragma once

#include "datastructures/SimulationConfiguration.h"
#include "datastructures/SimulationResult.h"
#include "ElasticScattering.h"

int sim_main(const InitParameters& init);
SimulationResult& RunSimulation(ElasticScattering& es, SimulationConfiguration& sp);
void PrintInfo(const SimulationConfiguration& sp, int count);
void LogResult(const SimulationConfiguration& sim_params, const SimulationResult& sr);
