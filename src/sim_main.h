#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "scattering/Grid.h"
#include "scattering/Simulation.h"

#include <string>

void sim_main(const InitParameters& init);

void RunSimulation(const SimulationConfiguration& sp, Simulation& es);

void UpdateMagneticField(ScatteringParameters& sp, double magnetic_field);
void UpdateTemperature(ScatteringParameters& sp, double temperature);
void CompleteSimulationParameters(ScatteringParameters& sp);

std::string GetAvailableDirectory(std::string base);
SimulationConfiguration ParseConfig(std::string file);
void PrintSimulationInfo(const SimulationConfiguration& cfg);