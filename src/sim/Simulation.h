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

	double GetTau(double temperature) const {
		return (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);
	}

	double GetDefaultMaxLifetime(double tau) const {
		return 15.0 * tau;
	}

	double GetSigmaIntegrandFactor(double tau) const {
		//	const double integrand_factor = sp->integrand_angle_area / ((values_per_quadrant - 1) * 3.0);
		double integral_factor = pow(1.0 / (3.0 * (ss.particles_per_row - 1)), 2);
		return integral_factor * SigmaFactor(tau) * 1e-8;
	}

	double AverageLifetime() const
	{
		if (ss.total_lifetimes == 0) return 0;

		int total_valid = 0;
		double total = 0;
		for (int i = 0; i < ss.total_lifetimes; i++) {
			total += raw_lifetimes[i];
			if (raw_lifetimes[i] > 0) total_valid++;
		}
		
		return total / (double)total_valid;
	}

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) = 0;
	virtual IterationResult DeriveTemperatureWithImages(const double temperature) const = 0;
	virtual Sigma			DeriveTemperature(const double temperature) const = 0;


	void InitSample(const Grid& grid, const Settings& s, const bool coherent)
	{
		is = grid.GetSettings();
		
		ss.region_size                = s.region_size;
		ss.region_extended_area       = s.region_extends;
		ss.distance_between_particles = ss.region_size / (double)(ss.particles_per_row - 1);
		ss.small_offset               = v2(ss.distance_between_particles * 0.01, ss.distance_between_particles * 0.005);

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

		//raw_lifetimes.resize(ss.particles_per_row * ss.values_per_row);
	}
};



class SimulationCPU : public Simulation {
public:
	SigmaResult ApplySigma(const double tau, const std::vector<double>& current_lifetimes) const;
	double IntegrateResult(const double tau, const std::vector<double>& sigma_lifetimes) const;

	SigmaResult ApplySigmaParticle(const double tau, const std::vector<double>& current_lifetimes) const;
	std::vector<double> IntegrateParticle(const std::vector<double>& current_lifetimes) const;
	double IntegrateSigma(const double tau, const std::vector<double>& particle_sigmas) const;

	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) override;
	virtual IterationResult DeriveTemperatureWithImages(const double temperature) const override;
	virtual Sigma			DeriveTemperature(const double temperature) const override;


	SimulationCPU(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant) {};
};


class SimulationCL : public Simulation {
	void PrepareKernels(const Settings& s, const size_t items_in_workgroup);

public:
	virtual void			ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics) override;
	virtual IterationResult DeriveTemperatureWithImages(const double temperature) const override;
	virtual Sigma			DeriveTemperature(const double temperature) const override;

	void UploadImpurities(const Grid& grid);

	SimulationCL();
	SimulationCL(int p_particles_per_row, int p_values_per_quadrant);
	// SimulationCL(bool use_gpu, bool show_info);
	~SimulationCL();
};

#endif // ELASTIC_SCATTERING_H
