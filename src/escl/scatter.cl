#include "src/escl/common.h"
#include "src/escl/util.h"

__kernel void lifetime(__global SimulationParameters *sp, __global double2 *impurities, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    if ((x < (row_size - 1)) && (y < (row_size - 1))) {
        DECLARE_POS

        lifetimes[y * row_size + x] = single_lifetime(pos, sp, impurities);
	}
    else { 
        lifetimes[y * row_size + x] = 0;
	}
}

__kernel void lifetime_phi(__global SimulationParameters *sp, __global double2 *impurities, __global double *integrand)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    
    DECLARE_POS

	integrand[y * row_size + x] = phi_lifetime(pos, sp, impurities);
}
