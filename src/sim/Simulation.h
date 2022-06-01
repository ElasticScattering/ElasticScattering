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

#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "Grid.h"

#include "es/settings.h"
#include "es/constants.h"

#include "src/ConfigSettings.h"
#include "src/SimulationResult.h"
#include "src/Metrics.h"

#include <vector>

struct PerformanceCounters
{
	LARGE_INTEGER lifetimeBegin, lifetimeEnd;
	LARGE_INTEGER temperaturesBegin, temperaturesEnd;
};

struct WorkSize
{
	size_t positions_per_row;
	size_t particles_per_position;
	int    total_particles;

	size_t particles_global[3];
	size_t particles_local [3];

	size_t sum_global;
	size_t sum_local;
	size_t summed_data_size;
};

// A Simulation operates on a generated impurity field.
// It can be run on the CPU, or on the GPU (OpenCL).
class Simulation {

protected:
	SimulationSettings ss;
	ImpuritySettings is;
	ParticleSettings ps;

	PerformanceCounters pc;

	LARGE_INTEGER clockFrequency;

	int GetIndex(int i, int j, int q, int p) const { return (j * ss.particles_per_row) + (i * ss.particles_per_position) + (q * ss.particles_per_quadrant) + p; };
	
	double SigmaFactor(double tau) const {
		double kf = M * ps.particle_speed / HBAR;
		double outside = (E * E * kf * kf) / (2.0 * PI * PI * M * C1);
		double wct = ps.angular_speed * tau;
		outside *= tau / (1.0 + wct * wct);

		return outside;
	}

	std::vector<double> GetTruncatedLifetimes(const double tau, const std::vector<double>& raw_lifetimes) const
	{
		const double default_max_lifetime = GetDefaultMaxLifetime(tau);
		std::vector<double> lifetimes(raw_lifetimes.size());
		for (int i = 0; i < lifetimes.size(); i++)
			lifetimes[i] = min(raw_lifetimes[i], default_max_lifetime);

		return lifetimes;
	}

	double GetSigmaIntegrandFactor(double tau) const { return pow(1.0 / (3.0 * (ss.positions_per_row - 1)), 2) * ss.phi_integrand_factor * SigmaFactor(tau) * 1e-8; }

	double GetTau(double temperature) const { return (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature); }
	double GetDefaultMaxLifetime(double tau) const { return 15.0 * tau; }

	double AverageLifetime(const std::vector<double> lifetimes) const
	{
		if (ss.total_particles == 0) return 0;

		int total_particles_started_outside_impurity = 0;
		double total = 0;
		for (int i = 0; i < ss.total_particles; i++) {
			total += lifetimes[i];
			
			if (lifetimes[i] > 0) total_particles_started_outside_impurity++;
		}
		
		return total / (double)ss.total_particles;
	}

	inline double GetElapsedTime(LARGE_INTEGER beginClock, LARGE_INTEGER endClock) const { return ((double)(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart); }

public:
	void InitSample(const Grid& grid, const Settings& s, const bool coherent)
	{
		is = grid.GetSettings();

		const double base_area  = s.alpha * 2.0;
		ss.integrand_angle_area = !coherent ? base_area : (HALF_PI - base_area);
		ss.phi_integrand_factor = ss.integrand_angle_area / ((ss.particles_per_quadrant - 1) * 3.0);
		ss.coherent_tau         = s.tau;

		ps.alpha          = s.alpha;
		ps.particle_speed = s.particle_speed;
		ps.is_clockwise   = s.is_clockwise ? 1 : 0;
		ps.is_coherent    = coherent ? 1 : 0;
		ps.phi_start      = (!coherent ? -s.alpha : s.alpha);
		ps.phi_step_size  = ss.integrand_angle_area / (ss.particles_per_quadrant - 1);
	}

	virtual std::vector<Sigma>			 ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) = 0;
	virtual std::vector<IterationResult> ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) = 0;

	Simulation(int p_positions_per_row, int p_particles_per_quadrant, const GridInformation& grid_info)
	{
		ss.positions_per_row      = p_positions_per_row;
		ss.particles_per_quadrant = p_particles_per_quadrant;
		ss.particles_per_position = ss.particles_per_quadrant * 4;
		ss.particles_per_row      = ss.positions_per_row * ss.particles_per_position;
		ss.total_particles        = ss.positions_per_row * ss.particles_per_row;
		ss.total_positions	      = ss.positions_per_row * ss.positions_per_row;

		ss.region_size                = grid_info.region_size;
		ss.region_extended_area       = grid_info.region_extends;
		ss.distance_between_positions = grid_info.region_size / (double)(ss.positions_per_row - 1);
		ss.small_offset               = v2(grid_info.cell_size * 0.01, grid_info.cell_size * 0.005);

		QueryPerformanceFrequency(&clockFrequency);
	}
};


class SimulationCPU : public Simulation {
	std::vector<double> ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics);
	Sigma			    DeriveTemperature(const double temperature, const std::vector<double> raw_lifetimes) const;
	IterationResult     DeriveTemperatureWithImages(const double temperature, const std::vector<double> raw_lifetimes) const;

	SigmaResult ApplySigma(const double tau, const std::vector<double>& current_lifetimes) const;
	double IntegrateResult(const double tau, const std::vector<double>& sigma_lifetimes) const;

	std::vector<double> IntegratePosition(const std::vector<double>& current_lifetimes) const;

public:
	virtual std::vector<Sigma>           ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;
	virtual std::vector<IterationResult> ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;

	SimulationCPU(int p_particles_per_row, int p_values_per_quadrant, const GridInformation& grid_info) : Simulation(p_particles_per_row, p_values_per_quadrant, grid_info) {};
};

class SimulationCL : public Simulation {
	WorkSize work_size;
	unsigned int last_grid_seed = 0;

	int coh_idx = 0;

	bool print_info = false;

	void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics);
	Sigma			DeriveTemperature(const double temperature) const;
	IterationResult DeriveTemperatureWithImages(const double temperature) const;

	void ApplyMaxLifetime(const double tau) const;
	void ApplySigmaComponent(const double tau) const;
	void ApplySimpsonWeights() const;
	Sigma SumBuffers(const double tau) const;
	
	void IntegrateToPositions() const;
	std::vector<double> ConvertToImage(std::vector<double> results) const;

public:
	virtual std::vector<Sigma>           ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;
	virtual std::vector<IterationResult> ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;

	SimulationCL(int p_particles_per_row, int p_values_per_quadrant, const GridInformation& grid_info, bool p_print_info);
	~SimulationCL();
};

#endif // ELASTIC_SCATTERING_H
