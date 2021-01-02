#include "src/sim/es/util.h"
#include "src/sim/cl/cl_macros.h"
#include "src/sim/es/constants.h"
#include "src/sim/es/settings.h"

kernel void 
apply_max_lifetime(constant double* raw_lifetimes, double default_max_lifetime, global double* lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);

	int values_per_particle = get_global_size(2);
	int particles_per_row   = get_global_size(0);
	int values_per_row      = particles_per_row * values_per_particle;

	int idx = j * values_per_row + i * values_per_particle + v;

	//int idx = GET_INDEX(i, j, v);
	double lt = raw_lifetimes[idx];
	lifetimes[idx] = min(lt , default_max_lifetime);
}


kernel void 
apply_sigma_component(constant double* lifetimes, 
					  constant SimulationSettings* ss, 
					  constant ParticleSettings* ps, 
					  const double tau, 
					  const int mode, 
					  global double* sigma_lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
	
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

	int idx    = GET_INDEX(i, j, v);
	double phi = GET_PHI(v);
	double f   = (mode == MODE_SIGMA_XX) ? cos(phi) : sin(phi); 
	sigma_lifetimes[idx] = GetSigma(lifetimes[idx], phi, tau, ss->signed_angular_speed) * f;
}


kernel void 
apply_simpson_weights(global double* sigma_lifetimes, constant SimulationSettings* ss, const int values_per_quadrant)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
    
	if (i >= ss->particles_per_row || j >= ss->particles_per_row) //@Todo: zo of met simpsonweight 0.
		return;

	int idx = GET_INDEX(i, j, v);

	sigma_lifetimes[idx] *= (SimpsonWeight2D(i, j, ss->particles_per_row) * SimpsonWeight(v % ss->values_per_quadrant, ss->values_per_quadrant));
}


kernel void 
integrate_to_particle(global double* values, constant SimulationSettings* ss, global double* particle_results)
{
	int i        = get_global_id(0);
    int j        = get_global_id(1);
	
	if (i >= ss->particles_per_row || j >= ss->particles_per_row) //@Todo: zo of met simpsonweight 0.
		return;
	
	int prt_idx = GET_PARTICLE_INDEX(i, j);
	double total = 0;
	for (int q = 0; q < 4; q++)
	{
		unsigned int base_idx = GET_INDEX(i, j, q); // @Optimize
	
		for (int p = 0; p < ss->values_per_quadrant; p++)
			total += values[base_idx + p] * SimpsonWeight(p, ss->values_per_quadrant);
	}
	
	// @Todo, integrand factor.
	particle_results[prt_idx] = total * ss->phi_integrand_factor;
}


kernel void 
sum(global double* A, global double* B, local double* local_sums)
{
	uint id = get_global_id(0);
	uint local_id = get_local_id(0);

	local_sums[local_id] = A[id] + A[id + get_global_size(0)]; // Global size should be half the length of A
	barrier(CLK_LOCAL_MEM_FENCE);
	for (int stride = get_local_size(0) / 2; stride > 1; stride /= 2) {
		if (local_id < stride) {
			local_sums[local_id] += local_sums[local_id + stride];
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (local_id == 0) {
		B[get_group_id(0)] = local_sums[0] + local_sums[1];
	}
}
