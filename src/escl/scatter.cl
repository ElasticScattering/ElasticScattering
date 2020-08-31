#include "src/escl/common.h"
#include "src/escl/util.h"

__kernel void lifetime(__global SimulationParameters *sp, __global double2 *impurities, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    if ((x < (row_size - 1)) && (y < (row_size - 1))) {
        //Remove 1 from row_size to have an inclusive range, another because the kernel work dimension is even, but the integral requires uneven dimensions.
        double2 pos = (double2)(x, y) * sp->region_size / (row_size-2); 
    
        double particle_lifetime;
        if (sp->angular_speed != 0) {
            bool clockwise = (sp->clockwise == 1);
            double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, clockwise, false);
            particle_lifetime = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, impurities);
	    } else {
            particle_lifetime = lifetime0(pos, sp, impurities);
        }

        lifetimes[y * row_size + x] = particle_lifetime;
	}
    else { 
        lifetimes[y * row_size + x] = 0;
	}
}

__kernel void sigma_xx(__global SimulationParameters *sp, __global double2 *impurities, __global double *integrand)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    
    double2 pos = (double2)(x, y) * (sp->region_size / (row_size-2));
    
    double angle_area = sp->alpha * 2.0;
    double step_size = angle_area / (sp->integrand_steps-1);
    double integral = 0;

    for (int j = 0; j < 4; j++)
    {
        double start = -sp->alpha + j * (PI * 0.5);
        double total = 0.0;
    
        bool is_even = true;
        for (int i = 0; i < sp->integrand_steps; i++)
        {
            sp->phi = start + i * step_size;

            double particle_lifetime;
            if (sp->angular_speed != 0) {
                bool clockwise = (sp->clockwise == 1);
                if (clockwise) sp->angular_speed *= -1;

                double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, clockwise, false);
                particle_lifetime = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, impurities);
			} else {
                particle_lifetime = lifetime0(pos, sp, impurities);
			}

	        double z = exp(-particle_lifetime / sp->tau);

            double r = cos(sp->phi) - cos(sp->phi + sp->angular_speed * particle_lifetime) * z;
	        r       += sp->angular_speed * sp->tau * sin(sp->phi + sp->angular_speed * particle_lifetime) * z;
	        r       -= sp->angular_speed * sp->tau * sin(sp->phi);
            r       *= sp->tau;
	        
            double rxx = r * cos(sp->phi);
            double rxy = r * sin(sp->phi);

            bool edge_item = (i == 0 || i == sp->integrand_steps-1);
            is_even = (i % 2) == 0;
            double w = 1.0;
    
            if (!edge_item) {
                w = is_even ? 2.0 : 4.0;
	        }

            total += rxx * w;
	    }

        integral += total * angle_area / ((sp->integrand_steps-1) * 3.0);
	}
    
	integrand[y * row_size + x] = integral;
}
