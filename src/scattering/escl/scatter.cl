#include "lifetime.h"
#include "constants.h"

__kernel void lifetime(__constant ScatteringParameters *sp, __global read_only double2 *impurities, __global read_only int *impurity_index, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    int limit = row_size-1;
    
    double result = 0;
    if ((x < limit) && (y < limit)) {
        //Remove 1 from row_size to have an inclusive range, another because the kernel work dimension is even, but the integral requires uneven dimensions.
        const double s = sp->region_size / (row_size - 2);
        const double2 pos = (double2)(x*s, y*s);

        if (sp->mode == MODE_DIR_LIFETIME) result = single_lifetime(pos, sp->phi, sp, impurities, impurity_index);
        else                               result = phi_lifetime   (pos,          sp, impurities, impurity_index);
	}

    lifetimes[y * row_size + x] = result;
}

/* ! Deprecated once impuritygrid works.
__kernel void scatter_sim(__constant ScatteringParameters *sp, __global double2 *impurities, __global double *xx, __global double *xy) 
{
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    uint row_size = get_global_size(0);
    uint limit = row_size-1;
    
    if ((x < limit) && (y < limit)) {
        const double2 pos = (double2)(x, y) * (sp->region_size / (row_size - 2));

        double2 result = SIM_phi_lifetime(pos, sp, impurities);

        uint idx = y * row_size + x;
        xx[idx] = result.x;
        xy[idx] = result.y;
    }
}
*/

__kernel void scatter_march(
            __constant ScatteringParameters *sp, 
            __global read_only double2 *impurities, __global read_only int *impurity_indices, 
            __global double *xx, __global double *xy)
{
    size_t x = get_global_id(0);
    uint y = get_global_id(1);
    uint row_size = get_global_size(0);
    uint limit = row_size-1;
    
    if ((x < limit) && (y < limit)) {
        const double2 pos = (double2)(x, y) * (sp->region_size / (row_size - 2));

        double2 result = march_phi_lifetime(pos, sp, impurities, impurity_indices);

        uint idx = y * row_size + x;
        xx[idx] = result.x;
        xy[idx] = result.y;
    }
}


__kernel void add_integral_weights_2d(__global double* A)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    int i = y * row_size + x;
    A[i] *= GetWeight2D(x, y, row_size-1);
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
