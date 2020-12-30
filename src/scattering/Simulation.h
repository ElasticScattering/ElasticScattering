#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "Grid.h"

#include "escl/ScatteringParameters.h"
#include "escl/constants.h"

#include "src/Settings.h"
#include "src/SimulationResult.h"
#include "src/Metrics.h"

#include <vector>

class Simulation {
protected:
	std::vector<double> raw_lifetimes;

	SimulationSettings ss;
	ImpuritySettings is;
	ParticleSettings ps;

	int GetIndex(int i, int j, int q, int p) const {
		return (j * ss.values_per_row) + (i * ss.values_per_particle) + (q * ss.values_per_quadrant) + p;
	};
	
	double SigmaFactor(double tau) const {
		double kf = M * ps.particle_speed / HBAR;
		double outside = (E * E * kf * kf) / (2.0 * PI * PI * M * C1);
		double wct = ps.angular_speed * tau;
		outside *= tau / (1.0 + wct * wct);

		return outside;
	}

	double AverageLifetime() const
	{
		if (ss.total_lifetimes == 0) return 0;

		double total = 0;
		for (int i = 0; i < ss.total_lifetimes; i++)
			total += raw_lifetimes[i];
		
		return total / (double)ss.total_lifetimes;
	}

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) = 0;
	virtual IterationResult DeriveTemperature(const double temperature) const = 0;

	void InitSample(const Grid& grid, const Settings& s, const bool coherent)
	{
		is = grid.GetSettings();
		
		ss.region_size                = s.region_size;
		ss.region_extended_area       = s.region_extends;
		ss.distance_between_particles = ss.region_size / (double)(ss.particles_per_row - 1);
		ss.small_offset               = v2(ss.distance_between_particles * 0.01, ss.distance_between_particles * 0.005);

		const double base_area = s.alpha * 2.0;
		ss.integrand_angle_area = !coherent ? base_area : (PI / 2.0 - base_area);
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

		//raw_lifetimes.resize(ss.particles_per_row * ss.values_per_row);
	}
};



class SimulationCPU : public Simulation {
public:
	SigmaResult ApplySigma(const double tau, const std::vector<double>& current_lifetimes) const;
	std::vector<double> IntegrateParticle(const std::vector<double>& current_lifetimes) const;
	double IntegrateSigma(const double tau, const std::vector<double>& particle_sigmas) const;

	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) override;
	virtual IterationResult DeriveTemperature(const double temperature) const override;

	SimulationCPU(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant) {};
};


class SimulationCL : public Simulation {
	void PrepareKernels(const Settings& ss, const size_t items_in_workgroup);

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) override;
	virtual IterationResult DeriveTemperature(const double temperature) const override;

	void UploadImpurities(const Grid& grid);

	SimulationCL();
	SimulationCL(int p_particles_per_row, int p_values_per_quadrant);;
	// SimulationCL(bool use_gpu, bool show_info);
	~SimulationCL();
};

#endif // ELASTIC_SCATTERING_H
