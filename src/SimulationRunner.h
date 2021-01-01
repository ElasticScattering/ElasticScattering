#pragma once

#include "SimulationConfiguration.h"
#include "scattering/Simulation.h"

#include <string>
#include <vector>

//Repeatedly starts simulations based on a config file, and logs the results.
class SimulationRunner {
	SimulationConfiguration cfg;

	LARGE_INTEGER clockFrequency;

	inline double GetElapsedTime(LARGE_INTEGER beginClock, LARGE_INTEGER endClock) const { return ((double)(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart); }

	std::string GetResultPath(int t_idx) const { return cfg.output_directory + "/T" + std::to_string(t_idx) + ".dat"; }
	std::string GetSamplePath(int sample_idx) const { return cfg.output_directory + "/Sample " + std::to_string(sample_idx); }

	std::string GetSampleResultsPath(int sample_idx) const { return cfg.output_directory + "/Sample " + std::to_string(sample_idx) + "/sample_results.dat"; }

	std::string GetImagePath(int t_idx, int m_idx, int sample_idx, bool coherent) const { 
		auto type = coherent ? "/Coherent/" : "/Incoherent/";
		return GetSamplePath(sample_idx) + type + "T" + std::to_string(t_idx) + " MF" + std::to_string(m_idx) + ".png";
	}

	//std::string GetMetricsPath(int sample_idx) const { return GetSamplePath(sample_idx) + "/Metrics.txt"; }

	std::string GetMetricsPath() const { return cfg.output_directory + "/Metrics.txt"; }

	void CreateOutputDirectory() const;
	void CreateSampleOutputDirectory(const int sample_index) const;

	SampleResult RunSample(Simulation& es, const Settings& settings, const int sample_index, const bool coherent, const Grid& grid);
	void FinishResults(const std::vector<SampleResult> sample_results_coh, const std::vector<SampleResult> sample_results_inc);

public:
	void Run(const InitParameters& init);
};
