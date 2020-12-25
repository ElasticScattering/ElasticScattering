#pragma once

#include "SimulationConfiguration.h"
#include "scattering/Simulation.h"

#include <string>

class SimulationRunner {
private:
	SimulationConfiguration cfg;

	LARGE_INTEGER beginClock, endClock, clockFrequency;

	inline double GetElapsedTime() { return ((double)(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart); }

	std::string GetSamplePath(int sample_idx) const;
	std::string GetImagePath(int t_idx, int m_idx, int sample_idx, bool incoherent) const;
	std::string GetMetricsPath(int m_idx, bool incoherent) const;
	std::string GetResultPath(int t_idx) const;

	void CompleteSimulationParameters(ScatteringParameters& sp);
	void UpdateMagneticField(ScatteringParameters& sp, double magnetic_field);
	void UpdateTemperature(ScatteringParameters& sp, double temperature);

	void ParseConfig(std::string file);
	std::string GetAvailableDirectory(std::string base);
	void CreateOutputDirectories() const;
	void PrintSimulationInfo() const;

	SampleResult RunSample(const int sample_index, Simulation& es, ScatteringParameters& sp, const Grid& grid);

public:
	void Run(const InitParameters& init);
};
