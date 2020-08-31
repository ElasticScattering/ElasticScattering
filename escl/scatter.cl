#include "escl/common.h"
#include "escl/util.h"

double lifetimeB(double max_lifetime, double2 pos, bool clockwise, __global SimulationParameters *sp, __global double2 *imps)
{
    double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = (double2)(cos(sp->phi), sin(sp->phi)) * sp->particle_speed;
    double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    double lifetime = max_lifetime;

    for (int i = 0; i < sp->impurity_count; i++) {
        double2 imp_pos = imps[i];

        double2 d = pos - imp_pos;

        if (CirclesCross(center, orbit_radius, imp_pos, sp->impurity_radius))
        {
            if (sp->impurity_radius_sq > dot(d,d))
            {
                lifetime = 0;
            }

            double t = GetFirstCrossTime(center, pos, imp_pos, orbit_radius, sp->impurity_radius, sp->angular_speed, clockwise);

            if (t < lifetime)
                lifetime = t;
		}
    }

    return lifetime;
}

double lifetime0(double max_lifetime, double2 pos, __global SimulationParameters *sp, __global double2 *imps)
{
    double2 unit = { cos(sp->phi), sin(sp->phi) };
    double2 vel = unit * sp->particle_speed;

    double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = sp->tau;
    
    for (int i = 0; i < sp->impurity_count; i++) {
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

__kernel void lifetime(__global SimulationParameters *sp, __global double2 *imps, __global double *lifetimes) 
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
            particle_lifetime = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, imps);
	    } else {
            particle_lifetime = lifetime0(sp->tau, pos, sp, imps);
        }

        lifetimes[y * row_size + x] = particle_lifetime;
	}
    else { 
        lifetimes[y * row_size + x] = 0;
	}
}

__kernel void sigma_xx(__global SimulationParameters *sp, __global double2 *imps, __global double *integrand)
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
                particle_lifetime = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, imps);
			} else {
                particle_lifetime = lifetime0(sp->tau, pos, sp, imps);
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
