#include "escl/common.h"

double smod(double a, double b)
{
    return a - b * floor(a / b);
}

inline double GetBoundTime(const double phi, const double alpha, const double w, const bool is_electron, const bool is_future)
{
    double remaining = smod(phi + alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2.0 * alpha - remaining;
    else if (is_future)            dphi = 2.0 * alpha - remaining;
    else                           dphi = remaining;

    return dphi / w;
}

double2 GetCyclotronOrbit(double2 p, double2 velocity, double radius, double vf, bool is_electron)
{
    double2 shift = (double2)(velocity.y, -velocity.x) * radius / vf;

    double2 center;
    if (is_electron)
    {
        center.x = p.x - shift.x;
        center.y = p.y - shift.y;
    }
    else { 
        center.x = p.x + shift.x;
        center.y = p.y + shift.y;
    }

    return center;
}

bool CirclesCross(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = q.x *q.x + q.y*q.y;
    if (dist_squared >= (r1 + r2)*(r1 + r2)) return false;
    if (dist_squared <= (r1 - r2)*(r1 - r2)) return false;

    return true;
}

double4 GetCrossPoints(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = dot(q,q);
    double dist = sqrt(dist_squared);
    double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    double ys = sqrt(r1 * r1 - xs * xs);

    double2 u = (double2)(p2 - p1) / dist;

    double4 points = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys,

        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

double GetPhi(double2 pos, double2 center, double radius)
{
    double p = (pos.x - center.x) / radius;
    if (p > 1.0) p = 1.0;
    if (p < -1.0) p = -1.0;
    double phi = acos(p);

    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

double GetCrossAngle(double p, double q, bool clockwise)
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

double GetFirstCrossTime(double2 center, double2 pos, double2 ip, double r, double ir, double w, double clockwise)
{
    double4 cross_points = GetCrossPoints(center, r, ip, ir);

    double2 p1 = { cross_points.x, cross_points.y };
    double2 p2 = { cross_points.z, cross_points.w };

    double phi0 = GetPhi(pos, center, r);
    double phi1 = GetPhi(p1, center, r);
    double phi2 = GetPhi(p2, center, r);

    double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;
    return min(t1, t2);
}

double lifetimeB(double max_lifetime, double2 pos, bool clockwise, SimulationParameters sp, __global double2 *imps)
{
    double orbit_radius = sp.particle_speed / sp.angular_speed;
    double2 vel = (double2)(cos(sp.phi), sin(sp.phi)) * sp.particle_speed;
    double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp.particle_speed, clockwise);

    double lifetime = max_lifetime;

    for (int i = 0; i < sp.impurity_count; i++) {
        double2 imp_pos = imps[i];

        double2 d = pos - imp_pos;

        if (CirclesCross(center, orbit_radius, imp_pos, sp.impurity_radius))
        {
            if (sp.impurity_radius_sq > dot(d,d))
            {
                lifetime = 0;
            }

            double t = GetFirstCrossTime(center, pos, imp_pos, orbit_radius, sp.impurity_radius, sp.angular_speed, clockwise);

            if (t < lifetime)
                lifetime = t;
		}
    }

    return lifetime;
}

double lifetime0(double max_lifetime, double2 pos, SimulationParameters sp, __global double2 *imps)
{
    double2 unit = { cos(sp.phi), sin(sp.phi) };
    double2 vel = unit * sp.particle_speed;

    double impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;

    double lifetime = sp.tau;
    
    for (int i = 0; i < sp.impurity_count; i++) {
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

    if ((x < (row_size - 1)) && (y < (row_size - 1))) {
        //Remove 1 from row_size to have an inclusive range, another because the kernel work dimension is even, but the integral requires uneven dimensions.
        double2 pos = (double2)(x, y) * sp.region_size / (row_size-1); 
    
        double particle_lifetime;
        if (sp.angular_speed != 0) {
            bool clockwise = (sp.clockwise == 1);
            double bound_time = GetBoundTime(sp.phi, sp.alpha, sp.angular_speed, clockwise, false);
            particle_lifetime = lifetimeB(min(sp.tau, bound_time), pos, clockwise, sp, imps);
	    } else {
            particle_lifetime = lifetime0(sp.tau, pos, sp, imps);
        }

        lifetimes[y * row_size + x] = particle_lifetime;
	}
}

__kernel void sigma_xx(SimulationParameters sp, __global double2 *imps, __global double *integrand)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    
    double2 pos = (double2)(x, y) * sp.region_size / (row_size-1);
    
    bool clockwise = (sp.clockwise == 1);
    if (clockwise) sp.angular_speed *= -1;

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

            double particle_lifetime;
            if (sp.angular_speed != 0) {
                double bound_time = GetBoundTime(phi, sp.alpha, sp.angular_speed, clockwise, false);
                particle_lifetime = lifetimeB(min(sp.tau, bound_time), pos, clockwise, sp, imps);
			}else {
                particle_lifetime = lifetime0(sp.tau, pos, sp, imps);
			}

	        double z = exp(-particle_lifetime / sp.tau);

            double r = cos(phi) - cos(phi + sp.angular_speed * particle_lifetime) * z;
	        r       += sp.angular_speed * sp.tau * sin(phi + sp.angular_speed * particle_lifetime) * z;
	        r       -= sp.angular_speed * sp.tau * sin(phi);
            r       *= sp.tau;
	        
            double rxx = r * cos(phi);
            double rxy = r * sin(phi);

            bool edge_item = (i == 0 || i == sp.integrand_steps-1);
            is_even = (i % 2) == 0;
            double w = 1.0;
    
            if (!edge_item) {
                w = is_even ? 2.0 : 4.0;
	        }
            

            total += rxx * w;
	    }

        integral += total * angle_area / ((sp.integrand_steps-1) * 3.0);
	}
    
	integrand[y * row_size + x] = integral;
}
