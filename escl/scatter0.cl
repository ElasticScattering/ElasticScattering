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

__kernel void lifetime(double region_size, double speed, double tau, double phi, int impurity_count, double imp_radius, __global double2 *imps, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    double2 pos = (double2)(region_size * x, region_size * y) / (row_size-1);

    lifetimes[y * row_size + x] = lifetime0(tau, pos, phi, speed, impurity_count, imp_radius, imps);
}
