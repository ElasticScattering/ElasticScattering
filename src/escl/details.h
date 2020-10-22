#pragma once
#include "constants.h"
#include "device_macros.h"

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
    const double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(r1 * r1 - xs * xs);

    const double2 u = (p2 - p1) / dist;

    double4 points = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys,

        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

inline double GetPhi(const double2 pos, const double2 center, const double radius)
{
    double p = (pos.x - center.x) / radius;

    p = (p > 1.0) ? 1.0 : p;
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