#include "src/scattering/escl/lifetime.h"
#include "src/scattering/escl/constants.h"

/*
__kernel void quadrant_lifetime(__constant ScatteringParameters* sp, __global read_only double2* impurities, __global read_only int* cell_indices)
{
	int x = get_global_id(0);
    int y = get_global_id(1);
	int q = get_global_id(2);

	int row_size = get_global_size(0);
    int limit = row_size-1;
	if (x >= limit || y >= limit) {
		return;
	}

	const double s = sp->region_size / (row_size - 2);
    const double2 pos = (double2)(x*s, y*s);

	double quadrant_start = q * (PI * 0.5);
	unsigned int index_base = y * limit + x * (sp->integrand_steps * 4) + (q*sp->integrand_steps);

	Particle particle;
	particle.starting_position = pos;
	
	int starting_cell = get_cell_index(pos, sp->impurity_spawn_range, sp->max_expected_impurities_in_cell);
	particle.cell_index = starting_cell;

	{
		// See if the particle starts inside an impurity.
		// If this is the case, we can skip all phi angles.
		int impurity_start = (p.cell_index - 1 < 0) ? 0 : cell_indices[p.cell_index - 1];
		int impurity_end = cell_indices[p.cell_index];

		bool starts_inside_impurity = false;
		for (int i = impurity_start; i < impurity_end; i++)
		{
			if (InsideImpurity(pos, impurities[i], sp->impurity_radius)) {
				starts_inside_impurity = true;
				break;
			}
		}

		if (starts_inside_impurity) {
			for (int i = 0; i < sp->integrand_steps; i++) {
				lifetimes[index_base+i] = 0;
			}

			return;
		}
	}

	// Trace the orbit for this quadrant.
	for (int i = 0; i < sp->integrand_steps; i++)
	{
		const double phi = sp->integrand_start_angle + quadrant_start + step * sp->integrand_step_size;
		const Orbit orbit = MakeOrbit(pos, phi, sp);
   
		particle.cell_index = starting_cell;
		particle.phi = phi;

		const double max_lifetime = min(sp->default_max_lifetime, orbit.bound_time);
		lifetimes[index_base+i] = TraceOrbit(&particle, &orbit, sp, impurities, cell_indices);
	}
}
*/
__kernel void add_simpson_weights_2d(__global double* A)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    int i = y * row_size + x;
    A[i] *= SimpsonWeight2D(x, y, row_size-1);
}

__kernel void sum(__global double* A, __global double* B, __local double* local_sums)
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

/*
#ifndef NO_WINDOW
__kernel void to_texture(__global double* lifetimes, int mode, double scale, __write_only image2d_t screen)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    float k = lifetimes[y * row_size + x] / scale;
    float4 c = (float4)(k, k, k, 1.0f);

    if (mode == MODE_SIGMA_XY) {
        if (k < 0.0) c = (float4)(0, 0, -k, 1.0f);
        else 		 c = (float4)(k, 0, 0, 1.0f);
    }

    write_imagef(screen, (int2)(x, y), c);
}
#endif // NO_WINDOW
*/