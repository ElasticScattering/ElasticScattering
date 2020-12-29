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

	double integrand_angle_area;
	double phi_integrand_factor;
	double signed_angular_speed; //?
	double coherent_tau;

	double region_size;
	double region_extended_area;

	v2 small_offset;

	int particles_per_row;
	int values_per_quadrant;
	int values_per_particle;
	int values_per_row;

	ImpuritySettings impurity_settings;
	ParticleSettings particle_settings;

	int GetIndex(int i, int j, int q, int p) const {
		return (j * values_per_row) + (i * values_per_particle) + (q * values_per_quadrant) + p;
	};
	
	double SigmaFactor(double tau, double region_size) const {
		double kf = M * particle_settings.particle_speed / HBAR;
		double outside = (E * E * kf * kf) / (2.0 * PI * PI * M * region_size * region_size * C1);
		double wct = particle_settings.angular_speed * tau;
		outside *= tau / (1.0 + wct * wct);

		return outside;
	}

	double AverageLifetime() const
	{
		const int count = raw_lifetimes.size();
		if (count == 0) return 0;

		double total = 0;
		for (int i = 0; i < count; i++) {
			total += raw_lifetimes[i];
		}
		return total / (double)count;
	}

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) = 0;
	virtual IterationResult DeriveTemperature(const double temperature) = 0;

	void InitSample(const Grid& grid, const Settings& ss, const bool coherent)
	{
		impurity_settings    = grid.GetSettings();
		region_size          = ss.region_size;
		region_extended_area = ss.region_extends;

		small_offset = v2(impurity_settings.cell_size * 0.01, impurity_settings.cell_size * 0.005);

		const double base_area = ss.alpha * 2.0;
		integrand_angle_area   = !coherent ? base_area : (PI / 2.0 - base_area);
		phi_integrand_factor   = integrand_angle_area / ((values_per_quadrant - 1) * 3.0);
		//voor volledige integratie: double factor = integrand_angle_area / (4 * (ss.integrand_steps - 1) * (limit * limit));

		particle_settings.particle_speed = ss.particle_speed;
		particle_settings.is_clockwise   = ss.is_clockwise ? 1 : 0;
		particle_settings.is_coherent    = coherent ? 1 : 0;
		particle_settings.phi_start      = (!coherent ? -ss.alpha : ss.alpha);
		particle_settings.phi_step_size  = integrand_angle_area / (values_per_quadrant - 1);

		coherent_tau = ss.tau;
	}

	Simulation(int p_particles_per_row, int p_values_per_quadrant)
	{
		particles_per_row   = p_particles_per_row;
		values_per_quadrant = p_values_per_quadrant;
		values_per_particle = values_per_quadrant * 4;
		values_per_row      = particles_per_row * values_per_particle;

		//raw_lifetimes.resize(particles_per_row * values_per_row);
	}
};



class SimulationCPU : public Simulation {
private:
	SigmaResult ApplySigma(const double tau, const std::vector<double>& current_lifetimes);
	std::vector<double> IntegrateParticle(const std::vector<double>& current_lifetimes);
	double IntegrateSigma(const double tau, const std::vector<double>& particle_sigmas);

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) override;
	virtual IterationResult DeriveTemperature(const double temperature) override;

	SimulationCPU(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant) {};
};


class SimulationCL : public Simulation {
	void PrepareKernels(const Settings& ss, const size_t items_in_workgroup);

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) override;
	virtual IterationResult DeriveTemperature(const double temperature) override;

	void UploadImpurities(const Grid& grid);

	SimulationCL();
	SimulationCL(int p_particles_per_row, int p_values_per_quadrant);;
	// SimulationCL(bool use_gpu, bool show_info);
	~SimulationCL();
};

#endif // ELASTIC_SCATTERING_H
