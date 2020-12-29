#include "Simulation.h"
#include "escl/lifetime.h"
#include "escl/util.h"

void SimulationCPU::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
	particle_settings.angular_speed = E * magnetic_field / M;
	signed_angular_speed            = particle_settings.is_clockwise ? -particle_settings.angular_speed : particle_settings.angular_speed;

	raw_lifetimes.clear();
	raw_lifetimes.resize(particles_per_row * values_per_row);

	for (int j = 0; j < particles_per_row; j++) {
		for (int i = 0; i < particles_per_row; i++) {
			const v2 pos = v2(i, j) * (region_size / (double)(particles_per_row - 1)) + small_offset;

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < values_per_quadrant; p++) {
					auto particle = CreateParticle(q, p, pos, &particle_settings);

					raw_lifetimes[GetIndex(i, j, q, p)] = TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), &metrics);
				}
			}
		}
	}

	metrics.avg_particle_lifetime = AverageLifetime();
}

IterationResult SimulationCPU::DeriveTemperature(const double temperature)
{
	double tau = (particle_settings.is_coherent) ? coherent_tau : HBAR / (KB * temperature);
	double default_max_lifetime = 15.0 * tau;

	IterationResult b;

	std::vector<double> new_lifetimes(raw_lifetimes.size());

	for (int i = 0; i < new_lifetimes.size(); i++)
		//new_lifetimes[i] = raw_lifetimes[i]; 
		new_lifetimes[i] = min(raw_lifetimes[i], default_max_lifetime);

	b.particle_lifetimes = IntegrateParticle(new_lifetimes);

	b.sigmas             = ApplySigma(tau, new_lifetimes);
	b.result.xx          = IntegrateSigma(tau, b.sigmas.xx_buffer);
	b.result.xy          = IntegrateSigma(tau, b.sigmas.xy_buffer);

	return b;
}

SigmaResult SimulationCPU::ApplySigma(const double tau, const std::vector<double>& current_lifetimes)
{
	std::vector<double> sigma_xx(particles_per_row * particles_per_row);
	std::vector<double> sigma_xy(particles_per_row * particles_per_row);

	double integral_total = 0;
	for (int j = 0; j < particles_per_row; j++) {
		for (int i = 0; i < particles_per_row; i++) {
			int particle_idx = j * particles_per_row + i;

			Sigma totals;

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < values_per_quadrant; p++) {
					double lt = current_lifetimes[GetIndex(i, j, q, p)];

					double phi = particle_settings.phi_start + q * (PI * 0.5) + p * particle_settings.phi_step_size;
					double sigma_base = GetSigma(lt, phi, tau, signed_angular_speed) * SimpsonWeight(p, values_per_quadrant);

					totals.xx += sigma_base * cos(phi);
					totals.xy += sigma_base * sin(phi);
				}
			}
			
			sigma_xx[particle_idx] = totals.xx * phi_integrand_factor;
			sigma_xy[particle_idx] = totals.xy * phi_integrand_factor;
		}
	}

	SigmaResult result;
	result.xx_buffer = sigma_xx;
	result.xy_buffer = sigma_xy;
	return result;
}

std::vector<double> SimulationCPU::IntegrateParticle(const std::vector<double>& current_lifetimes)
{
	std::vector<double> particle_lifetimes(particles_per_row * particles_per_row);

	for (int j = 0; j < particles_per_row; j++) {
		for (int i = 0; i < particles_per_row; i++) {
			double particle_total = 0;
			
			for (int q = 0; q < 4; q++) {
				for (size_t p = 0; p < values_per_quadrant; p++) {
					particle_total += current_lifetimes[GetIndex(i, j, q, p)] * SimpsonWeight(p, values_per_quadrant);
				}
			}

			int particle_idx = j * particles_per_row + i;
			particle_lifetimes[particle_idx] = particle_total * phi_integrand_factor;
		}
	}
	
	return particle_lifetimes;
}

double SimulationCPU::IntegrateSigma(const double tau, const std::vector<double>& particle_sigmas)
{
	double integral_total = 0;
	for (int j = 0; j < particles_per_row; j++) {
		double wy = SimpsonWeight(j, particles_per_row);
		for (int i = 0; i < particles_per_row; i++) {
			integral_total += particle_sigmas[j * particles_per_row + i] * (wy * SimpsonWeight(i, particles_per_row));
		}
	}
	
	double integral_factor = pow(1.0 / (3.0 * (particles_per_row - 1)), 2);
	double factor = integral_factor * SigmaFactor(tau);
	return integral_total * factor * 1e-8;
}
