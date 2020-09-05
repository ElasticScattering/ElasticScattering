#include "src/escl/common.h"

__kernel void lifetime(__global SimulationParameters *sp, __global double2 *impurities, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    if ((x < (row_size - 1)) && (y < (row_size - 1))) {
        //Remove 1 from row_size to have an inclusive range, another because the kernel work dimension is even, but the integral requires uneven dimensions.
        double s = sp->region_size / (row_size - 2);
        double2 pos = (double2)(x*s, y*s);

        double result;
        if (sp->mode == MODE_DIR_LIFETIME) result = single_lifetime(pos, sp, impurities);
        else                               result = phi_lifetime   (pos, sp, impurities);

        lifetimes[y * row_size + x] = result;
	}
    else { 
        lifetimes[y * row_size + x] = 0;
	}
}
