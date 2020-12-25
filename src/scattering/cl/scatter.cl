#include "src/scattering/escl/lifetime.h"
//#include "src/scattering/escl/constants.h"


#define GET_BASE_INDEX(i, j, q) j * ((sp->dim-1) * sp->values_per_particle) + (i * sp->values_per_particle) + (q * sp->integrand_steps)
#define GET_INDEX(i, j, q, p) GET_BASE_INDEX(i, j, q) + p
#define GET_PRT_INDEX(i, j) j * limit + i

#define GET_PHI(q, p) sp->integrand_start_angle + q * HALF_PI + p * sp->integrand_step_size


kernel void quadrant_lifetime(constant ScatteringParameters* sp, global double2* impurities, global int* imp_index, global double* lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int q = get_global_id(2);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || j >= limit)
		return;

	const double s = sp->region_size / (row_size - 2); 
    const double2 pos = ((double2)(i, j) * s) + (double2)(sp->cell_size * 0.01, sp->cell_size * 0.005);

	unsigned int base_idx = GET_BASE_INDEX(i, j, q);

	for (int p = 0; p < sp->integrand_steps; p++)
	{
		lifetimes[base_idx + p] = lifetime(q, p, pos, sp, impurities, imp_index);
	}
}

kernel void 
quadrant_apply_sigma_component(global double* lifetimes, constant ScatteringParameters* sp, int mode, global double* sigma_component)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int q = get_global_id(2);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || j >= limit)
		return;

	const double s = sp->region_size / (row_size - 2); 
    const double2 pos = ((double2)(i, j) * s) + (double2)(sp->cell_size * 0.01, sp->cell_size * 0.005);

	unsigned int base_idx = GET_BASE_INDEX(i, j, q);

	const double w = (sp->clockwise == 1) ? -sp->angular_speed : sp->angular_speed;

	for (int p = 0; p < sp->integrand_steps; p++)
	{
		double phi = GET_PHI(q, p);

		double f = (mode == MODE_SIGMA_XX) cos(p) : sin(p); 
		sigma_component[base_idx + p] = GetSigma(lifetimes[base_idx + p], p, sp->tau, w);
	}
}

kernel void 
integrate_to_particle(global double* values, constant ScatteringParameters* sp, global double* particle_results)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || y >= limit) {
		return;
	}

	double total = 0;
	for (int q = 0; q < 4; q++)
	{
		unsigned int base_idx = GET_BASE_INDEX(i, j, q);
	
		for (int p = 0; p < sp->integrand_steps; p++)
			total += values[base_idx + p] * SimpsonWeight(p, sp->integrand_steps);
	}
	
	const double integrand_factor = sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0);
	particle_results[GET_PRT_INDEX(i,j)] = totals.x * integrand_factor;
}

kernel void 
apply_simpson_weights_particles(global double* particle_values, constant ScatteringParameters* sp, global double* particle_results)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || y >= limit) {
		return;
	}

	const int idx = GET_PRT_INDEX(i,j);
	particle_results[idx] = particle_values[idx] * SimpsonWeight2D(x, y, row_size-1);
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
