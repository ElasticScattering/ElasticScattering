#include "ElasticScattering.h"
#include "escl/lifetime.h"
#include "escl/util.h"

std::vector<double> ElasticScatteringCPU::ComputeLifetimes(const ScatteringParameters& sp, const ImpurityIndex& grid)
{
	// GPU kernel works only with even work size.
	const int limit = sp.dim - 1;
	const int values_per_particle = sp.integrand_steps * 4;

	std::vector<double> lifetimes(limit * limit * values_per_particle);

	for (int j = 0; j < limit; j++) {
		printf("j: %d", j);
		
		for (int i = 0; i < limit; i++) {
			v2 pos = v2(i, j) * (sp.region_size / (double)(sp.dim - 2));

			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < sp.integrand_steps; p++) {
					lifetimes[j * limit + i * values_per_particle + (q * sp.integrand_steps) + p] = lifetime(q, p, pos, &sp, grid.GetImpurities(), grid.GetIndex());
				}
			}
		}
	}

	return lifetimes;
}

SigmaBuffer ElasticScatteringCPU::ComputeSigmas(const ScatteringParameters& sp, const std::vector<double>& lifetimes)
{
	const double w = (sp.is_clockwise == 1) ? -sp.angular_speed : sp.angular_speed;
	const int limit = sp.dim - 1;
	const int values_per_particle = sp.integrand_steps * 4;
	
	//@Todo, check of dit correct is.
	double integrand_factor = sp.integrand_angle_area / ((sp.integrand_steps - 1) * 3.0);

	std::vector<double> sigma_xx(limit * limit);
	std::vector<double> sigma_xy(limit * limit);

	double integral_total = 0;
	for (int j = 0; j < limit; j++) {
		for (int i = 0; i < limit; i++) {
			int particle_idx = j * limit + i;
			double particle_total_xx = 0;
			double particle_total_xy = 0;

			for (int q = 0; q < 4; q++) {
				int base_idx = j * limit + i * values_per_particle + (q * sp.integrand_steps);

				for (int p = 0; p < sp.integrand_steps; p++) {
					double lt = lifetimes[base_idx + p];

					double phi = sp.integrand_start_angle + q * (PI * 0.5) + p * sp.integrand_step_size;
					double sigma_base = GetSigma(lt, phi, sp.tau, w) * SimpsonWeight(p, sp.integrand_steps);

					particle_total_xx += sigma_base * cos(phi);
					particle_total_xy += sigma_base * sin(phi);
				}
			}
			
			sigma_xx[particle_idx] = particle_total_xx * integrand_factor;
			sigma_xy[particle_idx] = particle_total_xy * integrand_factor;
		}
	}

	SigmaBuffer buffer(sigma_xx, sigma_xy);
	return buffer;
}

std::vector<double> ElasticScatteringCPU::IntegrateParticle(const ScatteringParameters& sp, const std::vector<double>& lifetimes)
{
	const int limit = sp.dim - 1;
	const int values_per_particle = sp.integrand_steps * 4;
	
	std::vector<double> particle_lifetimes(lifetimes.size());

	for (int j = 0; j < limit; j++) {
		for (int i = 0; i < limit; i++) {
			int particle_idx = j * limit + i;
			double particle_total = 0;
			
			for (int q = 0; q < 4; q++) {
				int base_idx = j * limit + i * values_per_particle + (q * sp.integrand_steps);

				for (int p = 0; p < sp.integrand_steps; p++) {
					particle_total += lifetimes[base_idx + p] * SimpsonWeight(p, sp.integrand_steps);
				}
			}
			
			//@Todo, check of dit correct is.
			particle_lifetimes[particle_idx] = particle_total * sp.integrand_angle_area / ((sp.integrand_steps - 1) * 3.0);
		}
	}
	
	return particle_lifetimes;
}

SigmaResult ElasticScatteringCPU::IntegrateResult(const ScatteringParameters& sp, const std::vector<double>& lifetimes)
{
	// Apply weights for integration.
	const double w = (sp.is_clockwise == 1) ? -sp.angular_speed : sp.angular_speed;
	const int limit = sp.dim - 1;
	const int values_per_particle = sp.integrand_steps * 4;

	SigmaResult sr;

	double integral_total = 0;
	for (int j = 0; j < limit; j++) {
		double wy = SimpsonWeight(j, limit);
		for (int i = 0; i < limit; i++) {
			double wx = SimpsonWeight(i, limit);
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < sp.integrand_steps; p++) {
					double wp = SimpsonWeight(p, sp.integrand_steps);

					double lt = lifetimes[j * limit + i * values_per_particle + (q * sp.integrand_steps) + p];
					double phi = sp.integrand_start_angle + q * (PI * 0.5) + p * sp.integrand_step_size;

					double sigma_base = GetSigma(lt, phi, sp.tau, w) * (wy * wx * wp);
					sr.xx += sigma_base * cos(phi);
					sr.xy += sigma_base * cos(phi);
				}
			}
		}
	}

	double factor = sp.integrand_angle_area / (4 * (sp.integrand_steps - 1) * (limit * limit));
	factor *= SigmaFactor(sp);
	sr.xx *= factor;
	sr.xy *= factor;

	return sr;
}
