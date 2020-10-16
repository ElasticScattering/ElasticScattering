#include "src/escl/common.h"

__kernel void lifetime(__constant ScatteringParameters *sp, __global double2 *impurities, __global double *lifetimes) 
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

        if (sp->mode == MODE_DIR_LIFETIME) result = single_lifetime(pos, sp->phi, sp, impurities);
        else                               result = phi_lifetime   (pos,          sp, impurities);
	}

    lifetimes[y * row_size + x] = result;
}

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

__kernel void scatter_march(
            __constant ScatteringParameters *sp, 
            __global read_only double2 *impurities, 
            __global read_only uint *impurity_indices, 
            __global double *xx, 
            __global double *xy) 
{
    uint x = get_global_id(0);
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

#ifndef NO_WINDOW
__kernel void to_texture(__global double* lifetimes, int mode, double scale, __write_only image2d_t screen)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    float k = GetColor((float)(lifetimes[y * row_size + x]), scale, mode);
    float4 c = (float4)(k, k, k, 1.0f);

    if (mode == MODE_SIGMA_XY) {
        if (k < 0.0) c = (float4)(0, 0, -k, 1.0f);
        else 		 c = (float4)(k, 0, 0, 1.0f);
    }

    write_imagef(screen, (int2)(x, y), c);
}
#endif // NO_WINDOW
