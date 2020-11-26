#include "ElasticScattering.h"
#include "escl/lifetime.h"
#include "escl/util.h"


IterationResult ElasticScatteringCPU::ComputeIteration(const ScatteringParameters& sp, const ImpurityIndex& grid)
{
	auto lifetimes = ComputeLifetimes(sp, grid);

	IterationResult b;
	b.particle_lifetimes = IntegrateParticle(sp, lifetimes);
	b.sigmas = ComputeSigmas(sp, lifetimes);
	b.result = IntegrateResult(sp, lifetimes); //@Todo, gebruik sigma lifetimes hier?
	
	return b;
}

std::vector<double> ElasticScatteringCPU::ComputeLifetimes(const ScatteringParameters& sp, const ImpurityIndex& grid)
{
	// GPU kernel works only with even work size.
	const int limit = sp.dim - 1;

	std::vector<double> particle_lifetimes(limit * limit * sp.values_per_particle, 1);

	for (int j = 0; j < limit; j++) {
		for (int i = 0; i < limit; i++) {
			v2 pos = v2(i, j) * (sp.region_size / (double)(sp.dim - 2));

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < sp.integrand_steps; p++) {
					particle_lifetimes[j * limit + i * sp.values_per_particle + (q * sp.integrand_steps) + p] = lifetime(q, p, pos, &sp, grid.GetImpurities(), grid.GetIndex());
				}
			}
		}
	}

	return particle_lifetimes;
}

SigmaResult ElasticScatteringCPU::ComputeSigmas(const ScatteringParameters& sp, const std::vector<double>& particle_lifetimes)
{
	const double w = (sp.is_clockwise == 1) ? -sp.angular_speed : sp.angular_speed;
	const int limit = sp.dim - 1;
	
	//@Todo, check of dit correct is.
	const double integrand_factor = sp.integrand_angle_area / ((sp.integrand_steps - 1) * 3.0);

	std::vector<double> sigma_xx(limit * limit);
	std::vector<double> sigma_xy(limit * limit);

	double integral_total = 0;
	for (int j = 0; j < limit; j++) {
		for (int i = 0; i < limit; i++) {
			int particle_idx = j * limit + i;

			Sigma totals;

			for (int q = 0; q < 4; q++) {
				int base_idx = j * limit + i * sp.values_per_particle + (q * sp.integrand_steps);

				for (int p = 0; p < sp.integrand_steps; p++) {
					double lt = particle_lifetimes[base_idx + p];

					double phi = sp.integrand_start_angle + q * (PI * 0.5) + p * sp.integrand_step_size;
					double sigma_base = GetSigma(lt, phi, sp.tau, w) * SimpsonWeight(p, sp.integrand_steps);

					totals.xx += sigma_base * cos(phi);
					totals.xy += sigma_base * sin(phi);
				}
			}
			
			sigma_xx[particle_idx] = totals.xx * integrand_factor;
			sigma_xy[particle_idx] = totals.xy * integrand_factor;
		}
	}

	SigmaResult result;
	result.xx_buffer = sigma_xx;
	result.xy_buffer = sigma_xy;
	return result;
}

std::vector<double> ElasticScatteringCPU::IntegrateParticle(const ScatteringParameters& sp, const std::vector<double>& lifetimes)
{
	const int limit = sp.dim - 1;
	const double integrand_factor = sp.integrand_angle_area / ((sp.integrand_steps - 1) * 3.0);

	std::vector<double> particle_lifetimes(limit * limit);

	for (int j = 0; j < limit; j++) {
		for (int i = 0; i < limit; i++) {
			int particle_idx = j * limit + i;
			double particle_total = 0;
			
			for (int q = 0; q < 4; q++) {
				int base_idx = j * limit + i * sp.values_per_particle + (q * sp.integrand_steps);

				for (int p = 0; p < sp.integrand_steps; p++) {
					particle_total += lifetimes[base_idx + p] * SimpsonWeight(p, sp.integrand_steps);
				}
			}
			
			particle_lifetimes[particle_idx] = particle_total * integrand_factor;
		}
	}
	
	return particle_lifetimes;
}

Sigma ElasticScatteringCPU::IntegrateResult(const ScatteringParameters& sp, const std::vector<double>& particle_lifetimes)
{
	const double w = (sp.is_clockwise == 1) ? -sp.angular_speed : sp.angular_speed;
	const int limit = sp.dim - 1;

	Sigma sigma;

	double integral_total = 0;
	for (int j = 0; j < limit; j++) {
		double wy = SimpsonWeight(j, limit);
		for (int i = 0; i < limit; i++) {
			double wx = SimpsonWeight(i, limit);
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < sp.integrand_steps; p++) {
					double wp = SimpsonWeight(p, sp.integrand_steps);

					double lt = particle_lifetimes[j * limit + i * sp.values_per_particle + (q * sp.integrand_steps) + p];
					double phi = sp.integrand_start_angle + q * (PI * 0.5) + p * sp.integrand_step_size;

					double sigma_base = GetSigma(lt, phi, sp.tau, w) * (wy * wx * wp);
					sigma.xx += sigma_base * cos(phi);
					sigma.xy += sigma_base * cos(phi);
				}
			}
		}
	}

	double factor = sp.integrand_angle_area / (4 * (sp.integrand_steps - 1) * (limit * limit));
	factor *= SigmaFactor(sp);
	sigma.xx *= factor / 1e8;
	sigma.xy *= factor / 1e8;

	return sigma;
}
