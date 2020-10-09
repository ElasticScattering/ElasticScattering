#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "src/SimulationParameters.h"

#ifdef DEVICE_PROGRAM
    #include "src/escl/util.h"

    #define BUFFER_ARGS __constant SimulationParameters* sp, __global double2* impurities
#else
    #include "math.h"
    #include <windows.h>
    #include "src/escl/constants.h"
    #define BUFFER_ARGS SimulationParameters *sp, std::vector<v2> &impurities
    #define min(a, b) ((a) < (b)) ? (a) : (b)
#endif
    
inline bool IsEdge(int i, int dim) {
    return (i == 0) || (i == (dim - 1)); 
}

inline double GetWeight1D(int i, int dim) {
    double w = 1.0;
    if (!IsEdge(i, dim))
    {
        w *= ((i % 2) == 0) ? 2.0 : 4.0;
    }

    return w;
}

inline double GetWeight2D(int i, int j, int dim) {
    double w = 1.0;
    if (!IsEdge(i, dim))
    {
        w *= ((i % 2) == 0) ? 2.0 : 4.0;
    }
    if (!IsEdge(j, dim))
    {
        w *= ((j % 2) == 0) ? 2.0 : 4.0;
    }

    return w;
}

inline bool ShouldComputeSigma(int m) {
    return (m == MODE_SIMULATION || m == MODE_SIGMA_XX || m == MODE_SIGMA_XY);
}

inline double smod(const double a, const double b)
{
    return a - b * floor(a / b);
}

inline double GetBoundTime(const double phi, const double alpha, const double w, const bool is_incoherent, const bool is_diag_region, const bool is_electron, const bool is_future)
{
    if (!is_incoherent) return DBL_MAX;
    
    const double v = is_diag_region ? (phi + alpha - PI / 4.0) : (phi + alpha);
    const double remaining = smod(v, PI * 0.5);
    
    const double remaining2 = 2.0 * alpha - remaining;

    double dphi = ((!is_electron && is_future) || (is_electron && !is_future)) ? remaining : remaining2;
    dphi /= w;

    return dphi;
}

inline double2 GetCyclotronOrbit(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron)
{
    double2 shift = { velocity.y, -velocity.x };
    shift = shift * radius / vf;

    return is_electron ? (p - shift) : (p + shift);
}

inline bool CirclesCross(const double2 p1, const double r1, const double2 p2, const double r2)
{
    const double2 q = p1 - p2;
    const double dist_squared = q.x * q.x + q.y * q.y;
    
    const double r_add = r1 + r2;
    const double r_min = r1 - r2;

    return (dist_squared >= r_add * r_add || dist_squared <= r_min * r_min) ? false : true;
}

inline double4 GetCrossPoints(const double2 p1, const double r1, const double2 p2, const double r2)
{
    const double2 q = p1 - p2;

    const double dist_squared = dot(q, q);
    const double dist = sqrt(dist_squared);
    const double xs = (dist_squared + r1*r1 - r2*r2) / (2.0 * dist);
    const double ys = sqrt(r1*r1 - xs*xs);

    const double2 u = (p2 - p1) / dist;

    double4 points = {
        p1.x + u.x * xs +  u.y  * ys,
        p1.y + u.y * xs + -u.x  * ys,

        p1.x + u.x * xs + -u.y  * ys,
        p1.y + u.y * xs +  u.x  * ys
    };

    return points;
}

inline double GetPhi(const double2 pos, const double2 center, const double radius)
{
    double p = (pos.x - center.x) / radius;

    p = (p >  1.0) ?  1.0 : p;
    p = (p < -1.0) ? -1.0 : p;
    
    double phi = acos(p);
    phi = (pos.y < center.y) ? PI2 - phi : phi;

    return phi;
}

inline double GetCrossAngle(const double p, const double q, const bool clockwise)
{
    const double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

inline double GetFirstCrossTime(const double2 center, const double2 pos, const double2 ip, const double r, const double ir, const double w, const double clockwise)
{
    const double4 cross_points = GetCrossPoints(center, r, ip, ir);

    const double2 p1 = { cross_points.x, cross_points.y };
    const double2 p2 = { cross_points.z, cross_points.w };

    const double phi0 = GetPhi(pos, center, r);
    const double phi1 = GetPhi(p1, center, r);
    const double phi2 = GetPhi(p2, center, r);

    const double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    const double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;

    double a = min(2, 2);

    return min(t1, t2);
}

inline double lifetimeB(const double max_lifetime, const double2 pos, const double phi, const bool clockwise, BUFFER_ARGS)
{
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = { cos(phi), sin(phi) };
    vel = vel * sp->particle_speed;
    const double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);
    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = max_lifetime;

    for (int i = 0; i < sp->impurity_count; i++) {
        const double2 impurity = impurities[i];

        double2 d = pos - impurity;
        if (impurity_radius_sq > dot(d, d))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, orbit_radius, impurity, sp->impurity_radius))
        {
            double t = GetFirstCrossTime(center, pos, impurity, orbit_radius, sp->impurity_radius, sp->angular_speed, clockwise);

            lifetime = (t < lifetime) ? t : lifetime;
        }
    }

    return lifetime;
}

inline double lifetime0(const double2 pos, const double phi, BUFFER_ARGS)
{
    const double2 unit = { cos(phi), sin(phi) };
    const double2 vel = unit * sp->particle_speed;

    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = 15.0 * sp->tau;

    for (int i = 0; i < sp->impurity_count; i++) {
        const double2 imp_pos = impurities[i];
        const double inner = (imp_pos.x - pos.x) * unit.x + (imp_pos.y - pos.y) * unit.y;
        const double2 projected = pos + unit * inner;

        const double2 d = projected - imp_pos;
        const double diff = impurity_radius_sq - dot(d, d);
        if (diff < 0.0) {
            continue;
        }

        double L = sqrt(diff);

        double2 time_taken;
        
        if (fabs(vel.x) > (fabs(vel.y) * 1e-9)) {
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

////////////////////////////////////////////////////
// Main functions, for single or integrated lifetime
////////////////////////////////////////////////////

inline double single_lifetime(const double2 pos, const double phi, BUFFER_ARGS) {
    if (sp->angular_speed != 0) {
        const bool clockwise = (sp->is_clockwise == 1);
        const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, (sp->is_incoherent == 1), (sp->is_diag_regions == 1), clockwise, false);
 
        return lifetimeB(min(15.0 * sp->tau, bound_time), pos, phi, clockwise, sp, impurities);
    }
    else {
        return lifetime0(pos, phi, sp, impurities);
    }
}

inline double sigma_multiplier(double lt, double phi, double tau, double w)
{
    double z = exp(-lt / tau);

    double r = cos(phi) - cos(phi + w * lt) * z;

    r += w * tau * sin(phi + w * lt) * z;
    r -= w * tau * sin(phi);

    return r;
}

inline double phi_lifetime(const double2 pos, BUFFER_ARGS)
{
    const bool clockwise = (sp->is_clockwise == 1);
    const double w = (clockwise ? -sp->angular_speed : sp->angular_speed);

    bool diag_regions = (sp->is_diag_regions == 1);
    bool incoherent   = (sp->is_incoherent == 1); 
    
    const double incoherent_area = sp->alpha * 2.0;
    const double angle_area = incoherent ? incoherent_area : (PI / 2.0 - incoherent_area);
    
    const double step_size = angle_area / (sp->integrand_steps - 1);

    double base = (incoherent ? -sp->alpha : sp->alpha);
    if (diag_regions) {
        base += (incoherent ? (PI / 4.0) : -(PI / 4.0));
    }

    double integral = 0;
    for (int j = 0; j < 4; j++)
    {
        const double start = base + j * (PI * 0.5);
        double total = 0.0;

        for (int i = 0; i < sp->integrand_steps; i++)
        {
            const double phi = start + i * step_size;

            double result = single_lifetime(pos, phi, sp, impurities);

            if (ShouldComputeSigma(sp->mode)) {
                double multiplier = sigma_multiplier(result, phi, sp->tau, w);

                if (sp->mode != MODE_SIMULATION) {
                    double v;
                    if (sp->mode == MODE_SIGMA_XX) v = cos(phi);
                    else                           v = sin(phi);

                    result = multiplier * v;
                }
            }
      
            double w = 1.0;
            bool is_edge_value = !(i == 0 || i == (sp->integrand_steps - 1));

            if (is_edge_value)
            {
                w = ((i % 2) == 0) ? 2.0 : 4.0;
            }

            total += result * w;
        }

        integral += total * angle_area / ((sp->integrand_steps - 1) * 3.0);
    }

    return integral;
}

inline float GetColor(double v, double scale, int mode) {
    if      (mode == MODE_DIR_LIFETIME) v /= scale;
    else if (mode == MODE_PHI_LIFETIME) v /= (scale * 15.0);
    else if (mode == MODE_SIGMA_XX)     v /= 3.0;
    else if (mode == MODE_SIGMA_XY)     v /= 0.5;

    return (float)v;
}

#ifdef DEVICE_PROGRAM
__kernel void add_integral_weights_2d(__global double* A)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    int i = y * row_size + x;
    A[i] *= GetWeight2D(x, y, row_size-1);
}

__kernel void to_texture(__global double* lifetimes, int mode, double scale, __write_only image2d_t screen)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0);

    float k = GetColor((float)(lifetimes[y * row_size + x]), scale, mode);
    float4 c = (float4)(k, k, k, 1.0f);

    if (mode == MODE_SIGMA_XY) {
        if (k < 0.0) c = (float4)(0, 0, -k, 1.0f);
        else 		 c = (float4)(k, 0, 0, 1.0f);
    }

    write_imagef(screen, (int2)(x, y), c);
}

#endif

#endif // CL_COMMON_H
