#pragma once

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "scattering/Grid.h"
#include "scattering/Simulation.h"

#include <string>

void sim_main(const InitParameters& init);

void RunSimulation(const SimulationConfiguration& sp, Simulation& es);

SimulationConfiguration ParseConfig(const std::string file);
std::string GetAvailableDirectory(std::string base);
void CreateLog(const SimulationConfiguration& cfg, double temperature);
void LogResult(const std::string file_path, const DataRow row);
void FinishLog(const std::string file_path, const double time_elapsed);
void LogImages(const std::string file, const int dim, const double scale, const IterationResult iteration);
