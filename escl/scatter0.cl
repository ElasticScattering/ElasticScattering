#include "escl/common.h"



__kernel void lifetime(SimulationParameters sp, __global double2 *imps, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    if ((x < row_size - 1) && (y < row_size - 1)) {
        //Remove 1 from row_size to have an inclusive range, another because the kernel work dimension is even, but the integral requires uneven dimensions.
        double2 pos = (double2)(x, y) * sp.region_size / (row_size-2); 
        lifetimes[y * row_size + x] = lifetime0(sp.tau, pos, sp.phi, sp.particle_speed, sp.impurity_count, sp.impurity_radius, imps);
	}
}

__kernel void sigma_xx(SimulationParameters sp, __global double2 *imps, __global double *integrand)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    
    double2 pos = (double2)(x, y) * sp.region_size / (row_size-1);

    double angle_area = sp.alpha * 2.0;
    double step_size = angle_area / (sp.integrand_steps-1);
    double integral = 0;

    for (int j = 0; j < 4; j++)
    {
        double start = -sp.alpha + j * (PI * 0.5);
        double total = 0.0;
    
        bool is_even = true;
        for (int i = 0; i < sp.integrand_steps; i++)
        {
            double phi = start + i * step_size;

            double lt = lifetime0(sp.tau, pos, phi, sp.particle_speed, sp.impurity_count, sp.impurity_radius, imps);

	        double z = exp(-lt / sp.tau);

            double r = cos(phi) - cos(phi) * z;
	        r       += sp.angular_speed * sp.tau * sin(phi) * z;
	        r       -= sp.angular_speed * sp.tau * sin(phi);
            r       *= sp.tau;
	        
            double rxx = r * cos(phi);
            double rxy = r * sin(phi);

            bool edge_item = (i == 0 || i == sp.integrand_steps-1);
            
            double w = 1.0;
            if (!edge_item) {
                w = ((i % 2) == 0) ? 2.0 : 4.0;
	        }
            
            total += rxx * w;
	    }

        integral += total * angle_area / ((sp.integrand_steps-1) * 3.0);
	}
    
	integrand[y * row_size + x] = integral;
}
