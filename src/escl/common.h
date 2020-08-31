#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "common_structs.h"

#define PI   3.141592653589793238463
#define PI2  6.283185307179586
#define M0   9.109e-31
#define E    1.602e-19
#define HBAR 1.055e-34
#define C1   1.15e-9

#ifdef DEVICE_PROGRAM
    #define BUFFER_ARGS __global SimulationParameters* sp, __global double2* impurities
#else
    #include "math.h"

    #define BUFFER_ARGS SimulationParameters* sp, double2* impurities

    struct v4 {
        double x, y, z, w;
    };

    struct v2 {
        double x, y;

        v2(double p_x = 0, double p_y = 0)
            : x(p_x), y(p_y)
        {
        }

        v2& operator=(const v2& a)
        {
            x = a.x;
            y = a.y;
            return *this;
        }

        v2 operator+(const v2& a) const
        {
            return v2(a.x + x, a.y + y);
        }
        
        v2 operator-(const v2& a) const
        {
            return v2(x - a.x, y - a.y);
        }

        v2 operator/(double a) const
        {
            return v2(x / a, y / a);
        }

        v2 operator*(double a) const
        {
            return v2(x * a, y * a);
        }
    };

    #define double2 v2
    #define double4 v4

    inline double dot(double2 a, double2 b) { return a.x*b.x + b.y*b.y; }
#endif


inline double smod(double a, double b)
{
    return a - b * floor(a / b);
}

inline double GetBoundTime(double phi, double alpha, double w, bool is_electron, bool is_future)
{
    double remaining = smod(phi + alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2.0 * alpha - remaining;
    else if (is_future)            dphi = 2.0 * alpha - remaining;
    else                           dphi = remaining;

    return dphi / w;
}

inline double2 GetCyclotronOrbit(double2 p, double2 velocity, double radius, double vf, bool is_electron)
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

inline bool CirclesCross(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = q.x * q.x + q.y * q.y;
    if (dist_squared >= (r1 + r2) * (r1 + r2)) return false;
    if (dist_squared <= (r1 - r2) * (r1 - r2)) return false;

    return true;
}

inline double4 GetCrossPoints(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = dot(q, q);
    double dist = sqrt(dist_squared);
    double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    double ys = sqrt(r1 * r1 - xs * xs);

    double2 u = (p2 - p1) / dist;

    double4 points = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys,

        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

inline double GetPhi(double2 pos, double2 center, double radius)
{
    double p = (pos.x - center.x) / radius;
    if (p > 1.0) p = 1.0;
    if (p < -1.0) p = -1.0;
    double phi = acos(p);

    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

inline double GetCrossAngle(double p, double q, bool clockwise)
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

inline double GetFirstCrossTime(double2 center, double2 pos, double2 ip, double r, double ir, double w, double clockwise)
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


inline double lifetime0(double2 pos, BUFFER_ARGS)
{
    double2 unit = { cos(sp->phi), sin(sp->phi) };
    double2 vel = unit * sp->particle_speed;

    double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = sp->tau;

    for (int i = 0; i < sp->impurity_count; i++) {
        double2 imp_pos = impurities[i];
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

inline double lifetimeB(double max_lifetime, double2 pos, bool clockwise, BUFFER_ARGS)
{
    double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = (cos(sp->phi), sin(sp->phi)) * sp->particle_speed;
    double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    double lifetime = max_lifetime;

    for (int i = 0; i < sp->impurity_count; i++) {
        double2 imp_pos = impurities[i];

        double2 d = pos - imp_pos;

        if (CirclesCross(center, orbit_radius, imp_pos, sp->impurity_radius))
        {
            if (sp->impurity_radius_sq > dot(d, d))
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

#endif // CL_COMMON_H