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
	Metrics metrics;

	QueryPerformanceCounter(&pc.lifetimeBegin);
	auto lifetimes = ComputeLifetimes(magnetic_field, grid, metrics);
	QueryPerformanceCounter(&pc.lifetimeEnd);
	metrics.time_elapsed_lifetimes = GetElapsedTime(pc.lifetimeBegin, pc.lifetimeEnd);
	metrics.real_particles = ss.total_particles - metrics.particle_metrics.particles_inside_impurity;

	std::ofstream file;
	file.open("verification/cpu_lt.txt");
	file << std::setprecision(12);
	for (int j = 0; j < ss.positions_per_row; j++) {
		file << std::endl << "R" << j << std::endl;
		for (int i = 0; i < ss.positions_per_row; i++) {
			file << std::endl;
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.particles_per_quadrant; p++) {
					int idx = j * ss.positions_per_row * ss.particles_per_position + i * ss.particles_per_position + q * ss.particles_per_quadrant + p;

					file << lifetimes[idx] << std::endl;
				}
			}
		}
	}


	QueryPerformanceCounter(&pc.temperaturesBegin);
	std::vector<Sigma> results(temperatures.size());
	for (int i = 0; i < temperatures.size(); i++)
		results[i] = DeriveTemperature(temperatures[i], lifetimes);

	QueryPerformanceCounter(&pc.temperaturesEnd);
	metrics.time_elapsed_temperatures = GetElapsedTime(pc.temperaturesBegin, pc.temperaturesEnd);

	sample_metrics.iteration_metrics.push_back(metrics);

	return results;
}


std::vector<IterationResult> SimulationCPU::ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
	Metrics metrics;

	QueryPerformanceCounter(&pc.lifetimeBegin);
	auto lifetimes = ComputeLifetimes(magnetic_field, grid, metrics);
	QueryPerformanceCounter(&pc.lifetimeEnd);
	metrics.time_elapsed_lifetimes = GetElapsedTime(pc.lifetimeBegin, pc.lifetimeEnd);
	metrics.real_particles = ss.total_particles - metrics.particle_metrics.particles_inside_impurity;

	QueryPerformanceCounter(&pc.temperaturesBegin);
	std::vector<IterationResult> results(temperatures.size());
	for (int i = 0; i < temperatures.size(); i++)
		results[i] = DeriveTemperatureWithImages(temperatures[i], lifetimes);

	QueryPerformanceCounter(&pc.temperaturesEnd);
	metrics.time_elapsed_temperatures = GetElapsedTime(pc.temperaturesBegin, pc.temperaturesEnd);

	sample_metrics.iteration_metrics.push_back(metrics);

	return results;
}

///////////////////
///////////////////
///////////////////

std::vector<double> SimulationCPU::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
	std::vector<double> lifetimes(ss.total_particles);
	
	ps.angular_speed = E * magnetic_field / M;
	ss.signed_angular_speed = ps.is_clockwise ? -ps.angular_speed : ps.angular_speed;

	ParticleMetrics pm;

	for (int j = 0; j < ss.positions_per_row; j++) {
		for (int i = 0; i < ss.positions_per_row; i++) {
			const v2 position = v2(i, j) * ss.distance_between_positions + ss.small_offset;

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.particles_per_quadrant; p++) {
					auto particle = CreateParticle(q, p, position, &ps);

					int last_cells_passed = pm.cells_passed;
					int last_intersections_tested = pm.impurities_tested;

					lifetimes[GetIndex(i, j, q, p)] = TraceOrbit(&particle, &is, grid.GetImpurities(), grid.GetIndex(), pm);

					int new_cells_passed = pm.cells_passed - last_cells_passed;
					if (pm.max_cells_passed < new_cells_passed)
						pm.max_cells_passed = new_cells_passed;

					int new_intersections_tested = pm.impurities_tested - last_intersections_tested;
					if (pm.max_impurities_tested < new_intersections_tested)
						pm.max_impurities_tested = new_intersections_tested;
				}
			}
		}
	}

	metrics.particle_metrics = pm;
	metrics.avg_particle_lifetime = AverageLifetime(lifetimes);

	return lifetimes;
}

Sigma SimulationCPU::DeriveTemperature(const double temperature, const std::vector<double> raw_lifetimes) const
{
	double tau = GetTau(temperature);

	auto lifetimes = GetTruncatedLifetimes(tau, raw_lifetimes);
	auto sigmas = ApplySigma(tau, lifetimes);

	Sigma result;
	result.xx = IntegrateResult(tau, sigmas.xx_buffer);
	result.xy = IntegrateResult(tau, sigmas.xy_buffer);

	return result;
}

IterationResult SimulationCPU::DeriveTemperatureWithImages(const double temperature, const std::vector<double> raw_lifetimes) const
{
	double tau = GetTau(temperature);
	auto lifetimes = GetTruncatedLifetimes(tau, raw_lifetimes);
	
	IterationResult ir;
	ir.particle_lifetimes = IntegratePosition(lifetimes);
	
	auto sigmas		    = ApplySigma(tau, lifetimes);
	ir.sigmas.xx_buffer = IntegratePosition(sigmas.xx_buffer);
	ir.sigmas.xy_buffer = IntegratePosition(sigmas.xy_buffer);
	ir.result.xx        = IntegrateResult(tau, sigmas.xx_buffer);
	ir.result.xy        = IntegrateResult(tau, sigmas.xy_buffer);

	return ir;
}

//////////////////
//////////////////
//////////////////

SigmaResult SimulationCPU::ApplySigma(const double tau, const std::vector<double>& lifetimes) const
{
	std::vector<double> sigma_xx(ss.total_particles);
	std::vector<double> sigma_xy(ss.total_particles);

	double integral_total = 0;
	for (int j = 0; j < ss.positions_per_row; j++) {
		for (int i = 0; i < ss.positions_per_row; i++) {
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.particles_per_quadrant; p++) {
					int idx = GetIndex(i, j, q, p);

					double phi = ps.phi_start + q * HALF_PI + p * ps.phi_step_size;
					double sigma_base = GetSigma(lifetimes[idx], phi, tau, ss.signed_angular_speed);
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
	for (int j = 0; j < ss.positions_per_row; j++) {
		double wy = SimpsonWeight(j, ss.positions_per_row);

		for (int i = 0; i < ss.positions_per_row; i++) {
			double wp = wy * SimpsonWeight(i, ss.positions_per_row);

			for (int q = 0; q < 4; q++)
				for (int p = 0; p < ss.particles_per_quadrant; p++) 
					integral_total += sigma_lifetimes[GetIndex(i, j, q, p)] * (wp * SimpsonWeight(p, ss.particles_per_quadrant));
		}
	}

	return integral_total * GetSigmaIntegrandFactor(tau);
}

std::vector<double> SimulationCPU::IntegratePosition(const std::vector<double>& current_lifetimes) const
{
	std::vector<double> position_lifetimes(ss.total_positions);

	for (int j = 0; j < ss.positions_per_row; j++) {
		for (int i = 0; i < ss.positions_per_row; i++) {
			double position_total = 0;
			
			for (int q = 0; q < 4; q++) {
				for (size_t p = 0; p < ss.particles_per_quadrant; p++) {
					position_total += current_lifetimes[GetIndex(i, j, q, p)] * SimpsonWeight(p, ss.particles_per_quadrant);
				}
			}

			position_lifetimes[(j * ss.positions_per_row + i)] = position_total * ss.phi_integrand_factor;
		}
	}
	
	return position_lifetimes;
}
