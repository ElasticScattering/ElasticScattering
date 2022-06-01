/* -------------------------------------------------------------------------
	This code is part of ElasticScattering.
	
	Copyright(C) 2022 Elastic Scattering developers

	This program is free software : you can redistribute it and /or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#pragma once

#include "SimulationConfiguration.h"
#include "sim/Simulation.h"

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

	std::string GetMetricsPath() const { return cfg.output_directory + "/Metrics.txt"; }

	void CreateOutputDirectory() const;
	void CreateSampleOutputDirectory(const int sample_index) const;

	SimulationResult RunSimulation() const;
	SampleResult RunSample(Simulation& es, const Settings& settings, const int sample_index, const bool coherent, const Grid& grid) const;
	void FinishResults(const SimulationResult& sr) const;

public:
	void Run(const InitParameters& init);
};
