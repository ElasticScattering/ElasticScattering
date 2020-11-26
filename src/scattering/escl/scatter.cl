#include "src/scattering/escl/lifetime.h"
//#include "src/scattering/escl/constants.h"

kernel void quadrant_lifetimes(
	constant ScatteringParameters* sp,
	global read_only double2* impurities,
	global read_only int* imp_index,
	global write_only double* lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int q = get_global_id(2);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || y >= limit) {
		return;
	}

	const double s = sp->region_size / (row_size - 2);
    const double2 pos = (double2)(i, y) * s;

	unsigned int base_idx = y * limit + i * (sp->integrand_steps * 4) + (q*sp->integrand_steps);
	for (int p = 0; p < sp->integrand_steps; p++)
	{
		lifetimes[base_idx+p] = lifetime(q, p, pos, sp, impurities, imp_index);
	}
}

// ! Deprecated
kernel void quadrant_sigma(
	global read_only double* lifetimes, 
	constant ScatteringParameters* sp,
	global write_only double* sigma_xx, 
	global write_only double* sigma_xy)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || j >= limit) {
		return;
	}

	double2 total;
	for (int q = 0; q < 4; q++)
	{
		unsigned int base_idx = j * limit + i * (sp->integrand_steps * 4) + (q*sp->integrand_steps);
		
		for (int p = 0; p < sp->integrand_steps; p++)
		{
			double lt = lifetimes[base_idx+p];

			double phi = sp->integrand_start_angle + q * (PI * 0.5) + p * sp->integrand_step_size;
			double sigma_base = GetSigma(lt, phi, sp->tau, w) * SimpsonWeight(p, sp->integrand_steps);

			totals.x += sigma_base * cos(phi);
			totals.y += sigma_base * sin(phi);
		}
	}
	
	const int particle_idx = j * limit + i;
	const double integrand_factor = sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0);
	sigma_xx[particle_idx] = totals.x * integrand_factor;
	sigma_xy[particle_idx] = totals.y * integrand_factor;
}

kernel void quadrant_sigma_lifetimes(
	global read_only double* lifetimes, 
	constant ScatteringParameters* sp,
	global write_only double* sigma_xx, 
	global write_only double* sigma_xy)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int q = get_global_id(2);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || j >= limit) {
		return;
	}

	unsigned int base_idx = j * limit + i * (sp->integrand_steps * 4) + (q*sp->integrand_steps);
		
	for (int p = 0; p < sp->integrand_steps; p++)
	{
		double lt = lifetimes[base_idx+p];

		double phi = sp->integrand_start_angle + q * (PI * 0.5) + p * sp->integrand_step_size;
		double sigma_base = GetSigma(lt, phi, sp->tau, w);

		sigma_xx[particle_idx] = sigma_base * cos(phi);
		sigma_xy[particle_idx] = sigma_base * sin(phi);
	}
}

kernel void integrate_particles(
	global read_only double* lifetimes,
	constant ScatteringParameters* sp, 
	global write_only double* particle_lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);

	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || y >= limit) {
		return;
	}

	double total;
	for (int q = 0; q < 4; q++)
	{
		unsigned int base_idx = j * limit + i * (sp->integrand_steps * 4) + (q*sp->integrand_steps);
		
		for (int p = 0; p < sp->integrand_steps; p++)
		{
			total += lifetimes[base_idx+p] * SimpsonWeight(p, sp->integrand_steps);
		}
	}
	
	const double integrand_factor = sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0);
	particle_lifetimes[j * limit + i] = totals.x * integrand_factor;
}


kernel void add_simpson_weights_particle(global read_only double* A, global write_only double* B)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
	int q = get_global_id(2);
    int row_size = get_global_size(0);

    int i = y * row_size + x; 
    A[i] *= SimpsonWeight2D(x, y, row_size-1);
}


kernel void sum(global read_only double* A, global write_only double* B, local double* local_sums)
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