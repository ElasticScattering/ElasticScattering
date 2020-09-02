#include "src/escl/common.h"

__kernel void lifetime(__global SimulationParameters *sp, __global double2 *impurities, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    if ((x < (row_size - 1)) && (y < (row_size - 1))) {
        DECLARE_POS

        double result;
        if (sp->mode == MODE_DIR_LIFETIME) result = single_lifetime(pos, sp, impurities);
        else                               result = phi_lifetime   (pos, sp, impurities);

        lifetimes[y * row_size + x] = result;
	}
    else { 
        lifetimes[y * row_size + x] = 0;
	}
}
