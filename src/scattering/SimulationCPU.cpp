#include "Simulation.h"
#include "escl/lifetime.h"
#include "escl/util.h"

#include "src/Logger.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip> 

void SimulationCPU::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
	ps.angular_speed = E * magnetic_field / M;
	ss.signed_angular_speed = ps.is_clockwise ? -ps.angular_speed : ps.angular_speed;

	raw_lifetimes.clear();
	raw_lifetimes.resize(ss.total_lifetimes);

	for (int j = 0; j < ss.particles_per_row; j++) {
		for (int i = 0; i < ss.particles_per_row; i++) {
			const v2 pos = v2(i, j) * ss.distance_between_particles + ss.small_offset;

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.values_per_quadrant; p++) {
					auto particle = CreateParticle(q, p, pos, &ps);
					
					int last_cells_passed = metrics.cells_passed;
					int last_intersections_tested = metrics.impurities_tested;
					
					raw_lifetimes[GetIndex(i, j, q, p)] = TraceOrbit(&particle, &is, grid.GetImpurities(), grid.GetIndex(), &metrics);
					
					int new_cells_passed = metrics.cells_passed - last_cells_passed;
					if (metrics.max_cells_passed < new_cells_passed)
						metrics.max_cells_passed = new_cells_passed;
					
					int new_intersections_tested = metrics.impurities_tested - last_intersections_tested;
					if (metrics.max_impurities_tested < new_intersections_tested)
						metrics.max_impurities_tested = new_intersections_tested;
				}
			}
		}
	}
	
	/*	
	if (!ps.is_coherent) {
		std::ofstream file;
		file.open("ESLogs/mf " + std::to_string(magnetic_field) + ".txt");

		for (int i = 0; i < ss.total_lifetimes; i++)
		{
			file << std::scientific << std::setprecision(15) << raw_lifetimes[i] << std::endl;
		}
	}
	*/

	metrics.avg_particle_lifetime = AverageLifetime();
}

IterationResult SimulationCPU::DeriveTemperatureWithImages(const double temperature) const
{
	double tau = (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);
	double default_max_lifetime = 15.0 * tau;

	IterationResult ir;

	std::vector<double> new_lifetimes(raw_lifetimes.size());

	for (int i = 0; i < new_lifetimes.size(); i++)
		//new_lifetimes[i] = raw_lifetimes[i]; 
		new_lifetimes[i] = min(raw_lifetimes[i], default_max_lifetime);

	ir.particle_lifetimes = IntegrateParticle(new_lifetimes);
	
	/*
	auto sigmas		    = ApplySigma(tau, new_lifetimes);
	ir.sigmas.xx_buffer = IntegrateParticle(sigmas.xx_buffer);
	ir.sigmas.xy_buffer = IntegrateParticle(sigmas.xy_buffer);
	ir.result.xx        = IntegrateSigma(tau, ir.sigmas.xx_buffer);
	ir.result.xy        = IntegrateSigma(tau, ir.sigmas.xy_buffer);
	*/
	ir.sigmas = ApplySigmaParticle(tau, new_lifetimes);
	ir.result.xx = IntegrateSigma(tau, ir.sigmas.xx_buffer);
	ir.result.xy = IntegrateSigma(tau, ir.sigmas.xy_buffer);

	return ir;
}

Sigma SimulationCPU::DeriveTemperature(const double temperature) const
{
	double tau = (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);
	double default_max_lifetime = 15.0 * tau;

	std::vector<double> new_lifetimes(raw_lifetimes.size());

	for (int i = 0; i < new_lifetimes.size(); i++)
		//new_lifetimes[i] = raw_lifetimes[i]; 
		new_lifetimes[i] = min(raw_lifetimes[i], default_max_lifetime);

	auto sigmas = ApplySigma(tau, new_lifetimes);

	Sigma result;
	result.xx = IntegrateResult(tau, sigmas.xx_buffer);
	result.xy = IntegrateResult(tau, sigmas.xy_buffer);

	return result;
}

SigmaResult SimulationCPU::ApplySigma(const double tau, const std::vector<double>& current_lifetimes) const
{
	std::vector<double> sigma_xx(ss.total_lifetimes);
	std::vector<double> sigma_xy(ss.total_lifetimes);

	double integral_total = 0;
	for (int j = 0; j < ss.particles_per_row; j++) {
		for (int i = 0; i < ss.particles_per_row; i++) {
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.values_per_quadrant; p++) {
					int idx = GetIndex(i, j, q, p);
					double lt = current_lifetimes[idx];

					double phi = ps.phi_start + q * HALF_PI + p * ps.phi_step_size;
					double sigma_base = GetSigma(lt, phi, tau, ss.signed_angular_speed) * SimpsonWeight(p, ss.values_per_quadrant);
					sigma_xx[idx] = sigma_base * cos(phi);
					sigma_xy[idx] = sigma_base * sin(phi);
				}
			}
		}
	}

	SigmaResult result;
	result.xx_buffer = sigma_xx;
	result.xy_buffer = sigma_xy;

	return result;
}

double SimulationCPU::IntegrateResult(const double tau, const std::vector<double>& sigma_lifetimes) const
{
	double integral_total = 0;
	for (int j = 0; j < ss.particles_per_row; j++) {
		double wy = SimpsonWeight(j, ss.particles_per_row);

		for (int i = 0; i < ss.particles_per_row; i++) {
			double wp = wy * SimpsonWeight(i, ss.particles_per_row);

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.values_per_quadrant; p++) {
					integral_total += sigma_lifetimes[GetIndex(i, j, q, p)] * (wp * SimpsonWeight(i, ss.particles_per_row));
				}
			}
		}
	}

	/*
		stepsize = integrationwidth / (length - 2)
		return stepsize / 3

		//phirange = np.pi / 2 - 2 * alpha if iscoh else 2 * alpha
		const = integration_const(region, dimension)**2 * integration_const(phirange, nr_phi)
		const /= region**2
		integralxx *= const
		integralxy *= const
	*/

	double integral_factor = pow(1.0 / (3.0 * (ss.particles_per_row - 1)), 2);// *ss.phi_integrand_factor;
	double factor = integral_factor * SigmaFactor(tau);
	return integral_total * factor * 1e-8;
}

SigmaResult SimulationCPU::ApplySigmaParticle(const double tau, const std::vector<double>& current_lifetimes) const
{
	std::vector<double> sigma_xx(ss.total_particles);
	std::vector<double> sigma_xy(ss.total_particles);

	double integral_total = 0;
	for (int j = 0; j < ss.particles_per_row; j++) {
		for (int i = 0; i < ss.particles_per_row; i++) {
			int particle_idx = j * ss.particles_per_row + i;

			Sigma totals;

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.values_per_quadrant; p++) {
					double lt = current_lifetimes[GetIndex(i, j, q, p)];

					double phi = ps.phi_start + q * HALF_PI + p * ps.phi_step_size;
					double sigma_base = GetSigma(lt, phi, tau, ss.signed_angular_speed) * SimpsonWeight(p, ss.values_per_quadrant);

					totals.xx += sigma_base * cos(phi);
					totals.xy += sigma_base * sin(phi);
				}
			}
			
			sigma_xx[particle_idx] = totals.xx * ss.phi_integrand_factor;
			sigma_xy[particle_idx] = totals.xy * ss.phi_integrand_factor;
		}
	}

	SigmaResult result;
	result.xx_buffer = sigma_xx;
	result.xy_buffer = sigma_xy;
	return result;
}

std::vector<double> SimulationCPU::IntegrateParticle(const std::vector<double>& current_lifetimes) const
{
	std::vector<double> particle_lifetimes(ss.total_particles);

	for (int j = 0; j < ss.particles_per_row; j++) {
		for (int i = 0; i < ss.particles_per_row; i++) {
			double particle_total = 0;
			
			for (int q = 0; q < 4; q++) {
				for (size_t p = 0; p < ss.values_per_quadrant; p++) {
					particle_total += current_lifetimes[GetIndex(i, j, q, p)] * SimpsonWeight(p, ss.values_per_quadrant);
				}
			}

			int particle_idx = j * ss.particles_per_row + i;
			particle_lifetimes[particle_idx] = particle_total * ss.phi_integrand_factor;
		}
	}
	
	return particle_lifetimes;
}

double SimulationCPU::IntegrateSigma(const double tau, const std::vector<double>& particle_sigmas) const
{
	double integral_total = 0;
	for (int j = 0; j < ss.particles_per_row; j++) {
		double wy = SimpsonWeight(j, ss.particles_per_row);
		for (int i = 0; i < ss.particles_per_row; i++) {
			integral_total += particle_sigmas[j * ss.particles_per_row + i] * (wy * SimpsonWeight(i, ss.particles_per_row));
		}
	}
	
	double integral_factor = pow(1.0 / (3.0 * (ss.particles_per_row - 1)), 2);
	double factor = integral_factor * SigmaFactor(tau);
	return integral_total * factor * 1e-8;
}
