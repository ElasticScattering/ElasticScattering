#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "common_structs.h"

#define PI   3.141592653589793238463
#define PI2  6.283185307179586
#define M0   9.109e-31
#define E    1.602e-19
#define HBAR 1.055e-34
#define C    1.15e-9

double smod(double a, double b)
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

    double dist_squared = q.x * q.x + q.y * q.y;
    if (dist_squared >= (r1 + r2) * (r1 + r2)) return false;
    if (dist_squared <= (r1 - r2) * (r1 - r2)) return false;

    return true;
}

double4 GetCrossPoints(double2 p1, double r1, double2 p2, double r2)
{
    double2 q = p1 - p2;

    double dist_squared = dot(q, q);
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

#endif // CL_COMMON_H