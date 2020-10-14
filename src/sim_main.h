#pragma once

#include "datastructures/SimulationConfiguration.h"
#include "datastructures/SimulationResult.h"
#include "ElasticScattering.h"

int sim_main(const InitParameters& init);
void ComputeIteration(SimulationElasticScattering& es, SimulationConfiguration& sp, SimulationResult& sr);
void PrintInfo(const SimulationConfiguration& sp, int count);
