#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "Grid.h"

#include "es/settings.h"
#include "es/constants.h"

#include "src/ConfigSettings.h"
#include "src/SimulationResult.h"
#include "src/Metrics.h"

#include <vector>

class Simulation {

protected:
	SimulationSettings ss;
	ImpuritySettings is;
	ParticleSettings ps;

	LARGE_INTEGER clockFrequency;

	int GetIndex(int i, int j, int q, int p) const { return (j * ss.values_per_row) + (i * ss.values_per_particle) + (q * ss.values_per_quadrant) + p; };
	
	double SigmaFactor(double tau) const {
		double kf = M * ps.particle_speed / HBAR;
		double outside = (E * E * kf * kf) / (2.0 * PI * PI * M * C1);
		double wct = ps.angular_speed * tau;
		outside *= tau / (1.0 + wct * wct);

		return outside;
	}

	double GetSigmaIntegrandFactor(double tau) const { return pow(1.0 / (3.0 * (ss.particles_per_row - 1)), 2) * ss.phi_integrand_factor * SigmaFactor(tau) * 1e-8; }

	double GetTau(double temperature) const { return (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature); }
	double GetDefaultMaxLifetime(double tau) const { return 15.0 * tau; }

	double AverageLifetime(const std::vector<double> lifetimes) const
	{
		if (ss.total_lifetimes == 0) return 0;

		int total_valid = 0;
		double total = 0;
		for (int i = 0; i < ss.total_lifetimes; i++) {
			total += lifetimes[i];
			if (lifetimes[i] > 0) total_valid++;
		}
		
		return total / (double)total_valid;
	}

	inline double GetElapsedTime(LARGE_INTEGER beginClock, LARGE_INTEGER endClock) const { return ((double)(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart); }

public:
	virtual std::vector<Sigma> ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) = 0;
	virtual std::vector<IterationResult> ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) = 0;

	void InitSample(const Grid& grid, const Settings& s, const bool coherent)
	{
		is = grid.GetSettings();
		
		ss.region_size                = s.region_size;
		ss.region_extended_area       = s.region_extends;
		ss.distance_between_particles = ss.region_size / (double)(ss.particles_per_row - 1);
		ss.small_offset = v2(is.cell_size * 0.01, is.cell_size * 0.005);

		const double base_area = s.alpha * 2.0;
		ss.integrand_angle_area = !coherent ? base_area : (HALF_PI - base_area);
		ss.phi_integrand_factor = ss.integrand_angle_area / ((ss.values_per_quadrant - 1) * 3.0);
		ss.coherent_tau         = s.tau;

		ps.alpha		  = s.alpha;
		ps.particle_speed = s.particle_speed;
		ps.is_clockwise   = s.is_clockwise ? 1 : 0;
		ps.is_coherent    = coherent ? 1 : 0;
		ps.phi_start      = (!coherent ? -s.alpha : s.alpha);
		ps.phi_step_size  = ss.integrand_angle_area / (ss.values_per_quadrant - 1);
	}

	Simulation(int p_particles_per_row, int p_values_per_quadrant)
	{
		ss.particles_per_row   = p_particles_per_row;
		ss.values_per_quadrant = p_values_per_quadrant;
		ss.values_per_particle = ss.values_per_quadrant * 4;
		ss.values_per_row      = ss.particles_per_row * ss.values_per_particle;
		ss.total_lifetimes     = ss.particles_per_row * ss.values_per_row;
		ss.total_particles	   = ss.particles_per_row * ss.particles_per_row;

		QueryPerformanceFrequency(&clockFrequency);

		//raw_lifetimes.resize(ss.particles_per_row * ss.values_per_row);
	}
};


class SimulationCPU : public Simulation {
	std::vector<double> ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics);
	Sigma			    DeriveTemperature(const double temperature, const std::vector<double> raw_lifetimes) const;
	IterationResult     DeriveTemperatureWithImages(const double temperature, const std::vector<double> raw_lifetimes) const;

	SigmaResult ApplySigma(const double tau, const std::vector<double>& current_lifetimes) const;
	double IntegrateResult(const double tau, const std::vector<double>& sigma_lifetimes) const;
	std::vector<double> IntegrateParticle(const std::vector<double>& current_lifetimes) const;

public:
	virtual std::vector<Sigma>           ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;
	virtual std::vector<IterationResult> ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;

	SimulationCPU(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant) {};
};


class SimulationCL : public Simulation {
	void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics);
	Sigma			DeriveTemperature(const double temperature) const;
	IterationResult DeriveTemperatureWithImages(const double temperature) const;

	void PrepareKernels(const Settings& s, const size_t items_in_workgroup);
public:
	virtual std::vector<Sigma>           ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;
	virtual std::vector<IterationResult> ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics) override;

	void UploadImpurities(const Grid& grid);

	SimulationCL(int p_particles_per_row, int p_values_per_quadrant);
	~SimulationCL();
};

#endif // ELASTIC_SCATTERING_H
