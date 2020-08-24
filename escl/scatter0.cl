#include "escl/common.h"

double lifetime0(double tau, double2 pos, double phi, double speed, int impurity_count, double imp_radius, __global double2 *imps)
{
    double impurity_radius_sq = imp_radius * imp_radius;
    double lifetime = tau;
      
    double2 unit = { cos(phi), sin(phi) };
    double2 vel = unit * speed;
    
    for (int i = 0; i < impurity_count; i++) {
        double2 imp_pos = imps[i];
        double inner = (imp_pos.x - pos.x) * unit.x + (imp_pos.y - pos.y) * unit.y;
        double2 projected = pos + unit * inner;
        
        double2 d = projected - imp_pos;
        double diff = impurity_radius_sq - dot(d, d);
        if (diff < 0.0) {
            continue;
        }

        double L = sqrt(diff);

        double2 time_taken;
        if (vel.x != 0) {
            time_taken.x = -((projected.x - L * unit.x) - pos.x) / vel.x;
            time_taken.y = -((projected.x + L * unit.x) - pos.x) / vel.x;
        }
        else {
            time_taken.x = -((projected.y - L * unit.y) - pos.y) / vel.y;
            time_taken.y = -((projected.y + L * unit.y) - pos.y) / vel.y;
        }

        if ((time_taken.x * time_taken.y) < 0) {
            lifetime = 0;
            break;        
		}

        if (time_taken.x > 0 && time_taken.x < lifetime) {
            lifetime = time_taken.x;
        }
        if (time_taken.y > 0 && time_taken.y < lifetime) {
            lifetime = time_taken.y;
        }
    }

    return lifetime;
}

__kernel void lifetime(SimulationParameters sp, __global double2 *imps, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    double2 pos = (double2)(x, y) * sp.region_size / (row_size-1);

    lifetimes[y * row_size + x] = lifetime0(sp.tau, pos, sp.phi, sp.particle_speed, sp.impurity_count, sp.impurity_radius, imps);
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
