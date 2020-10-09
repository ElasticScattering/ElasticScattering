#include "src/escl/common.h"

__kernel void lifetime(__constant SimulationParameters *sp, __global double2 *impurities, __global double *lifetimes) 
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


__kernel void scatter_sim(__constant SimulationParameters *sp, __global double2 *impurities, __global double *xx, __global double *xy) 
{
    uint x = get_global_id(0);
    uint y = get_global_id(1);
    uint row_size = get_global_size(0);
    uint limit = row_size-1;
    
    double result = 0;
    if ((x < limit) && (y < limit)) {
        const double s = sp->region_size / (row_size - 2);
        const double2 pos = (double2)(x*s, y*s);

        result = phi_lifetime(pos, sp, impurities);
	}

    uint idx = y * row_size + x;
    xx[idx] = result * cos(phi);
    xy[idx] = result * sin(phi);
}
