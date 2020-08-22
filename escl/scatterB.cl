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

double lifetime0(double max_lifetime, double2 pos, double phi, double speed, double angular_speed, bool clockwise, int impurity_count, double imp_radius,  __global double2 *imps)
{
    double lifetime = max_lifetime;

    double radius = speed / angular_speed;
    double2 vel = (double2)(cos(phi), sin(phi)) * speed;
    double2 center = GetCyclotronOrbit(pos, vel, radius, speed, clockwise);

    double impurity_radius_sq = imp_radius * imp_radius;

    for (int i = 0; i < impurity_count; i++) {
        double2 imp_pos = imps[i];

        double2 d = pos - imp_pos;
        if (impurity_radius_sq > dot(d,d))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, radius, imp_pos, imp_radius))
        {
            double t = GetFirstCrossTime(center, pos, imp_pos, radius, imp_radius, angular_speed, clockwise);
            if (t < lifetime)
                lifetime = t;
		}
    }

    return lifetime;
}
__kernel void lifetime(SimulationParameters sp, __global double2 *imps, __global double *lifetimes) 
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    
    double2 pos = (double2)(sp.region_size * x, sp.region_size * y) / (row_size-1);
    
    bool clockwise = false;
    double bound_time = GetBoundTime(sp.phi, sp.alpha, sp.angular_speed, clockwise, false);
    
    lifetimes[y * row_size + x] = lifetime0(min(sp.tau, bound_time), pos, sp.phi, sp.particle_speed, sp.angular_speed, clockwise, sp.impurity_count, sp.impurity_radius, imps);
}

__kernel void sigma_xx(SimulationParameters sp, __global double2 *imps, __global double *integrand)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);
    
    double2 pos = (double2)(sp.region_size * x, sp.region_size * y) / (row_size-1);
    bool clockwise = false;
    if (clockwise) sp.angular_speed *= -1;

    int steps = 49;
    double angle_area = sp.alpha * 2.0;
    double step_size = angle_area / (steps-1);
    double integral = 0;

    for (int j = 0; j < 4; j++)
    {
        double start = -sp.alpha + j * (PI * 0.5);
        double total = 0.0;
    
        bool is_even = true;
        for (int i = 0; i < steps; i++)
        {
            double phi = start + i * step_size;

            double bound_time = GetBoundTime(phi, sp.alpha, sp.angular_speed, clockwise, false);
    
            double lt = lifetime0(min(sp.tau, bound_time), pos, phi, sp.particle_speed, sp.angular_speed, clockwise, sp.impurity_count, sp.impurity_radius, imps);

	        double z = exp(-lt / sp.tau);

            double r = cos(phi) - cos(phi + sp.angular_speed * lt) * z;
	        r       += sp.angular_speed * sp.tau * sin(phi + sp.angular_speed * lt) * z;
	        r       -= sp.angular_speed * sp.tau * sin(phi);
            r       *= sp.tau;
	        
            double rxx = r * cos(phi);
            double rxy = r * sin(phi);

            bool edge_item = (i == 0 || i == steps-1);
            //is_even = !is_even; // faster
            is_even = (i % 2) == 0;
            double w;
    
            if (edge_item) {
                w = 1.0;
	        } else {
                w = is_even ? 2.0 : 4.0;
	        }
            

            total += rxx * w;
	    }

        integral += total * angle_area / ((steps-1) * 3.0);
	}
    
	integrand[y * row_size + x] = integral;
}

