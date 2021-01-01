#include "src/sim/es/util.h"
#include "shared_macros.h"

kernel void 
apply_max_lifetime(constant double* raw_lifetimes, double default_max_lifetime, global double* lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);

	int idx = GET_INDEX(i, j, v);
	double lt = raw_lifetimes[idx];
	lifetimes[idx] = min(lt , default_max_lifetime);
}

kernel void 
apply_sigma_component(constant double* lifetimes, constant SimulationSettings* ss, int mode, global double* sigma_component)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
	
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

	unsigned int idx = GET_INDEX(i, j, v);
	const int q = (int)(v / ss->integrand_steps);
	const int p =       v % ss->integrand_steps;

	//double tau = (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);

	double f = (mode == MODE_SIGMA_XX) cos(p) : sin(p); 
	sigma_component[idx] = GetSigma(lifetimes[idx], GET_PHI(q, p), tau, ss->signed_angular_speed);
}

kernel void 
apply_simpson_weights(global double* sigma_lifetimes, double values_per_quadrant, global double* B)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
	int limit = get_global_size(0) - 1;
    
	if (i >= limit || y >= limit) //@Todo: zo of met simpsonweight 0.
		return;

	unsigned int idx = GET_INDEX(i, j, v);
	B[idx] = sigma_lifetimes[idx] * (SimpsonWeight2D(i, j, limit) * SimpsonWeight(v % values_per_quadrant, values_per_quadrant));
}

kernel void 
integrate_to_particle(global double* values, double values_per_quadrant, double integrand_factor, global double* particle_results)
{
	int i        = get_global_id(0);
    int j        = get_global_id(1);
	int row_size = get_global_size(0) - 1;
	int limit    = row_size - 1;
	
	if (i >= limit || y >= limit) //@Todo: zo of met simpsonweight 0.
		return;

	double total = 0;
	for (int q = 0; q < 4; q++)
	{
		unsigned int base_idx = GET_INDEX(i, j, q); // @Optimize
	
		for (int p = 0; p < values_per_quadrant; p++)
			total += values[base_idx + p] * SimpsonWeight(p, values_per_quadrant);
	}
	
	particle_results[j * row_size + i] = total;
}

kernel void 
sum(global read_only double* A, global write_only double* B, local double* local_sums)
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
