#pragma once
#include "details.h"

#include "constants.h"

#ifndef DEVICE_PROGRAM
    #include <cmath>
    #include <cfloat>
#endif


double smod(const double a, const double b)
{
    return a - b * floor(a / b);
}

double GetBoundTime(const double phi, const double alpha, const double w, const bool is_incoherent, const bool is_diag_region, const bool is_electron, const bool is_future)
{
    if (!is_incoherent) return DBL_MAX;

    const double v = is_diag_region ? (phi + alpha - PI / 4.0) : (phi + alpha);
    const double remaining = smod(v, PI * 0.5);

    const double remaining2 = 2.0 * alpha - remaining;

    double dphi = ((!is_electron && is_future) || (is_electron && !is_future)) ? remaining : remaining2;
    dphi /= w;

    return dphi;
}

double2 GetCyclotronOrbitCenter(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron)
{
    double2 shift = { velocity.y, -velocity.x };
    shift = shift * radius / vf;

    return is_electron ? (p - shift) : (p + shift);
}

bool CirclesCross(const Orbit c1, const double2 p2, const double r2)
{
    const double2 q = c1.position - p2;
    const double dist_squared = q.x * q.x + q.y * q.y;

    const double r_add = c1.radius + r2;
    const double r_min = c1.radius - r2;

    return (dist_squared >= r_add * r_add || dist_squared <= r_min * r_min) ? false : true;
}

double4 GetCrossPoints(const Orbit c1, const double2 p2, const double r2)
{
    const double2 q = c1.position - p2;

    const double dist_squared = dot(q, q);
    const double dist = sqrt(dist_squared);
    const double xs = (dist_squared + c1.radius_squared - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(c1.radius_squared - xs * xs);

    const double2 u = (p2 - c1.position) / dist;

    double4 points = {
        c1.position.x + u.x * xs + u.y * ys,
        c1.position.y + u.y * xs + -u.x * ys,

        c1.position.x + u.x * xs + -u.y * ys,
        c1.position.y + u.y * xs + u.x * ys
    };

    return points;
}

double GetPhi(const double2 pos, const Orbit c1)
{
    double p = (pos.x - c1.position.x) / c1.radius;

    p = (p > 1.0) ? 1.0 : p;
    p = (p < -1.0) ? -1.0 : p;

    double phi = acos(p);
    phi = (pos.y < c1.position.y) ? PI2 - phi : phi;

    return phi;
}

double GetCrossAngle(const double p, const double q, const bool clockwise)
{
    const double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

double GetFirstCrossTime(const double2 pos, const Orbit c1, const double2 ip, const double ir, const double w)
{
    const double4 cross_points = GetCrossPoints(c1, ip, ir);

    const double2 p1 = { cross_points.x, cross_points.y };
    const double2 p2 = { cross_points.z, cross_points.w };

    const double phi0 = GetPhi(pos, c1);
    const double phi1 = GetPhi(p1, c1);
    const double phi2 = GetPhi(p2, c1);

    const double t1 = GetCrossAngle(phi0, phi1, c1.clockwise) / w;
    const double t2 = GetCrossAngle(phi0, phi2, c1.clockwise) / w;

    return min(t1, t2);
}
