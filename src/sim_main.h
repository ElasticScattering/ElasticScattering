#pragma once

int sim_main(const InitParameters& init);
void ComputeIteration(SimulationElasticScattering& es, SimulationParameters& sp, SimulationResult& sr);
void PrintInfo(const SimulationParameters& sp, int count);
