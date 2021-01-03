#include "Simulation.h"
#include "es/lifetime.h"
#include "es/util.h"

#include "src/Logger.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <Windows.h>

std::vector<Sigma> SimulationCPU::ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
	LARGE_INTEGER beginLifetimesClock, endLifetimesClock;

	Metrics metrics;

	QueryPerformanceCounter(&beginLifetimesClock);
	auto lifetimes = ComputeLifetimes(magnetic_field, grid, metrics);
	QueryPerformanceCounter(&endLifetimesClock);
	metrics.time_elapsed_lifetimes = GetElapsedTime(beginLifetimesClock, endLifetimesClock);
	
	metrics.real_lifetimes = ss.total_lifetimes - metrics.particles_inside_impurity;
	sample_metrics.iteration_metrics.push_back(metrics);

	std::vector<Sigma> results(temperatures.size());
	for (int i = 0; i < temperatures.size(); i++)
		results[i] = DeriveTemperature(temperatures[i], lifetimes);

	return results;
}

std::vector<IterationResult> SimulationCPU::ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
	LARGE_INTEGER beginLifetimesClock, endLifetimesClock;
	Metrics metrics;

	QueryPerformanceCounter(&beginLifetimesClock);
	auto lifetimes = ComputeLifetimes(magnetic_field, grid, metrics);
	QueryPerformanceCounter(&endLifetimesClock);
	metrics.time_elapsed_lifetimes = GetElapsedTime(beginLifetimesClock, endLifetimesClock);

	metrics.real_lifetimes = ss.total_lifetimes - metrics.particles_inside_impurity;
	sample_metrics.iteration_metrics.push_back(metrics);

	std::vector<IterationResult> results(temperatures.size());
	for (int i = 0; i < temperatures.size(); i++)
		results[i] = DeriveTemperatureWithImages(temperatures[i], lifetimes);

	return results;
}

///////////////////
///////////////////
///////////////////

std::vector<double> SimulationCPU::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
	std::vector<double> lifetimes(ss.total_lifetimes);
	
	ps.angular_speed = E * magnetic_field / M;
	ss.signed_angular_speed = ps.is_clockwise ? -ps.angular_speed : ps.angular_speed;

	for (int j = 0; j < ss.particles_per_row; j++) {
		for (int i = 0; i < ss.particles_per_row; i++) {
			const v2 pos = v2(i, j) * ss.distance_between_particles + ss.small_offset;

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.values_per_quadrant; p++) {
					auto particle = CreateParticle(q, p, pos, &ps);

					int last_cells_passed = metrics.cells_passed;
					int last_intersections_tested = metrics.impurities_tested;

					lifetimes[GetIndex(i, j, q, p)] = TraceOrbit(&particle, &is, grid.GetImpurities(), grid.GetIndex(), &metrics);

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

	metrics.avg_particle_lifetime = AverageLifetime(lifetimes);

	return lifetimes;
}

Sigma SimulationCPU::DeriveTemperature(const double temperature, const std::vector<double> raw_lifetimes) const
{
	double tau = GetTau(temperature);
	double default_max_lifetime = GetDefaultMaxLifetime(tau);

	std::vector<double> new_lifetimes(raw_lifetimes.size());
	for (int i = 0; i < new_lifetimes.size(); i++)
		new_lifetimes[i] = min(raw_lifetimes[i], default_max_lifetime);

	auto sigmas = ApplySigma(tau, new_lifetimes);

	Sigma result;
	result.xx = IntegrateResult(tau, sigmas.xx_buffer);
	result.xy = IntegrateResult(tau, sigmas.xy_buffer);

	return result;
}

IterationResult SimulationCPU::DeriveTemperatureWithImages(const double temperature, const std::vector<double> raw_lifetimes) const
{
	double tau = GetTau(temperature);
	double default_max_lifetime = GetDefaultMaxLifetime(tau);

	std::vector<double> new_lifetimes(raw_lifetimes.size());
	for (int i = 0; i < new_lifetimes.size(); i++)
		new_lifetimes[i] = min(raw_lifetimes[i], default_max_lifetime);
	
	IterationResult ir;
	ir.particle_lifetimes = IntegrateParticle(new_lifetimes);
	
	auto sigmas		    = ApplySigma(tau, new_lifetimes);
	ir.sigmas.xx_buffer = IntegrateParticle(sigmas.xx_buffer);
	ir.sigmas.xy_buffer = IntegrateParticle(sigmas.xy_buffer);
	ir.result.xx        = IntegrateResult(tau, sigmas.xx_buffer);
	ir.result.xy        = IntegrateResult(tau, sigmas.xy_buffer);

	return ir;
}

//////////////////
//////////////////
//////////////////

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
					double sigma_base = GetSigma(lt, phi, tau, ss.signed_angular_speed);
					sigma_xx[idx] = sigma_base *cos(phi);
					sigma_xy[idx] = sigma_base *sin(phi);
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
					integral_total += sigma_lifetimes[GetIndex(i, j, q, p)] * (wp * SimpsonWeight(p, ss.values_per_quadrant));
				}
			}
		}
	}

	return integral_total * GetSigmaIntegrandFactor(tau);
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

			particle_lifetimes[(j * ss.particles_per_row + i)] = particle_total * ss.phi_integrand_factor;
		}
	}
	
	return particle_lifetimes;
}
