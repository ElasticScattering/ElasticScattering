#pragma once

#pragma once

#include "SimulationConfiguration.h"
#include "scattering/Simulation.h"

#include <string>


class SimulationRunner {
private:
	SimulationConfiguration cfg;

	std::string GetImagePathBase(int t_idx, int m_idx) const;
	std::string GetImagePath(int t_idx, int m_idx, int s_idx, bool incoherent) const;
	std::string GetMetricsPath(int m_idx) const;
	std::string GetResultPath(int t_idx) const;

	void CompleteSimulationParameters(ScatteringParameters& sp);
	void UpdateMagneticField(ScatteringParameters& sp, double magnetic_field);
	void UpdateTemperature(ScatteringParameters& sp, double temperature);

	void ParseConfig(std::string file);
	std::string GetAvailableDirectory(std::string base);
	void CreateOutputDirectories() const;
	void PrintSimulationInfo() const;

	void RunIteration(const int magnetic_index, Simulation& es);

public:
	void Run(const InitParameters& init);
};