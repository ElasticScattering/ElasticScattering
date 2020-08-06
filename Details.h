#ifndef DETAILS_H
#define DETAILS_H

#include "Constants.h"
#include <assert.h>

inline double smod(double a, double b)
{
    return a - b * floor(a / b);
}

inline double GetBoundTime(const double phi, const double alpha, const double w, const bool is_electron, const bool is_future)
{
    double remaining = smod(phi + alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2 * alpha - remaining;
    else if (is_future)            dphi = 2 * alpha - remaining;
    else                           dphi = remaining;

    return dphi / w;
}

inline v2 GetCyclotronOrbit(const v2 p, const v2 velocity, const double radius, const double vf, const bool is_electron)
{
    v2 shift = { radius * velocity.y / vf, -radius * velocity.x / vf };

    v2 center;
    if (is_electron) center = { p.x - shift.x, p.y - shift.y };
    else             center = { p.x + shift.x, p.y + shift.y };

    return center;
}

inline bool CirclesCross(const v2 p1, const double r1, const v2 p2, const double r2)
{
    double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    if (dist_squared >= pow(r1 + r2, 2)) return false;
    if (dist_squared <= pow(r1 - r2, 2)) return false;

    return true;
}

inline std::pair<v2, v2> GetCrossPoints(const v2 p1, const double r1, const v2 p2, const double r2)
{
    const double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    const double dist = sqrt(dist_squared);
    const double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(r1 * r1 - xs * xs);

    v2 u = { (p2.x - p1.x) / dist, (p2.y - p1.y) / dist };

    std::pair<v2, v2> points;
    points.first = {
        p1.x + u.x * xs + u.y * ys,
        p1.y + u.y * xs + -u.x * ys
    };

    points.second = {
        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs + u.x * ys
    };

    return points;
}

inline double GetPhi(const v2 pos, const v2 center, const double radius)
{
    double p = (pos.x - center.x) / radius;
    assert(abs(p) < 1.0001);
    p = max(min(p, 1), -1);
    double phi = acos(p);

    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

inline double GetCrossAngle(const double p, const double q, const bool clockwise)
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

inline double GetFirstCrossTime(const v2 center, const v2 pos, const v2 ip, const double r, const double ir, const double w, const double clockwise)
{
    const auto cross_points = GetCrossPoints(center, r, ip, ir);

    const double phi0 = GetPhi(pos, center, r);
    const double phi1 = GetPhi(cross_points.first, center, r);
    const double phi2 = GetPhi(cross_points.second, center, r);

    const double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    const double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;
    return MIN(t1, t2);
}

#endif // DETAILS_H