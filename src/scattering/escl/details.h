#pragma once

#include "constants.h"
#include "device_macros.h"

#ifndef DEVICE_PROGRAM
    #include <cmath>
    #include <cfloat>
#endif

typedef struct Orbit {
	double2 center;
	double radius;
	double radius_squared;
	
	double bound_time;
	double bound_phi;

	bool clockwise;

	Orbit(double2 _center, double _radius, bool _clockwise, double _bound_time, double _bound_phi) {
		center = _center;
		radius = _radius;
		clockwise = _clockwise;
		bound_time = _bound_time;
		bound_phi = _bound_phi;

		radius_squared = radius * radius;
	}
} Orbit;

typedef struct Particle {
	double phi;
	double2 starting_position;
    int2 starting_cell;
} Particle;


ESCL_INLINE double smod(const double a, const double b)
{
    return a - b * floor(a / b);
}

// @Optimize
ESCL_INLINE bool AngleInRange(const double phi, const double2 phi_range, bool clockwise)
{
    if (abs(phi_range.x - phi_range.y) < 1e-10) return true;

    double low = (clockwise) ? phi_range.y : phi_range.x;
    double high = (clockwise) ? phi_range.x : phi_range.y;

    high = (high < low) ? (high + PI2) : high;

    return (low < phi && phi < high ) || (low < (phi + PI2) && (phi + PI2) < high);
}

ESCL_INLINE bool CirclesCross(const Orbit* orbit, const double2 p2, const double r2)
{
    const double2 q = orbit->center - p2;
    const double dist_squared = q.x * q.x + q.y * q.y;

    const double r_add = orbit->radius + r2;
    const double r_min = orbit->radius - r2;

    return (dist_squared < r_add * r_add) && (dist_squared > r_min * r_min);
}

ESCL_INLINE double4 GetCrossPoints(const Orbit* orbit, const double2 p2, const double r2)
{
    const double2 q = orbit->center - p2;

    const double dist_squared = dot(q, q);
    const double dist = sqrt(dist_squared);
    const double xs = (dist_squared + orbit->radius_squared - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(orbit->radius_squared - xs * xs);

    const double2 u = (p2 - orbit->center) / dist;

    double4 points = {
        orbit->center.x + u.x * xs + u.y * ys,
        orbit->center.y + u.y * xs + -u.x * ys,

        orbit->center.x + u.x * xs + -u.y * ys,
        orbit->center.y + u.y * xs + u.x * ys
    };

    return points;
}

ESCL_INLINE double GetPhi(const double2 pos, const Orbit* orbit)
{
    double p = (pos.x - orbit->center.x) / orbit->radius;

    p = (p > 1.0) ? 1.0 : p;
    p = (p < -1.0) ? -1.0 : p;

    double phi = acos(p);
    phi = (pos.y < orbit->center.y) ? PI2 - phi : phi;

    return phi;
}

ESCL_INLINE double GetPositionAngle(const double angle_velocity, const bool clockwise)
{
    double offset = clockwise ? -PI / 2 : PI / 2;

    return fmod(angle_velocity + offset, PI2);
}

ESCL_INLINE double GetCrossAngle(const double p, const double q, const bool clockwise)
{
    const double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

/* Returns the first intersection with the impurity circle, or INF if it lies in an invalid range.
 * 
 * Orbit hits this impurity, but it can still lie in a subsection of the cell that is strictly later on 
 * the orbit's trajectory than parts of other cells, and thus it should be ignored until we have found 
 * no intersection in those cells.
 */
ESCL_INLINE double GetFirstCrossTime(const Orbit* orbit, const double2 pos, const double2 ip, const double ir, const double w, const double2 valid_range)
{
    const double4 cross_points = GetCrossPoints(orbit, ip, ir);

    const double2 p1 = { cross_points.x, cross_points.y };
    const double2 p2 = { cross_points.z, cross_points.w };

    const double phi0 = GetPhi(pos, orbit);
    const double phi1 = GetPhi(p1, orbit);
    const double phi2 = GetPhi(p2, orbit);

    double traversal_time = INF;

    if (AngleInRange(phi1, valid_range, orbit->clockwise)) {
        const double t = GetCrossAngle(phi0, phi1, orbit->clockwise) / w;
        traversal_time = min(traversal_time, t);
    }

    if (AngleInRange(phi2, valid_range, orbit->clockwise)) {
        const double t = GetCrossAngle(phi0, phi2, orbit->clockwise) / w;
        traversal_time = min(traversal_time, t);
    }

    return traversal_time;
}

ESCL_INLINE double GetBoundTime(const double phi, const double alpha, const double w, const bool is_incoherent, const bool is_diag_region, const bool is_electron, const bool is_future)
{
    if (!is_incoherent) return INF;

    const double v = is_diag_region ? (phi + alpha - PI / 4.0) : (phi + alpha);
    const double remaining = smod(v, PI * 0.5);

    // Oud
    //double dphi = ((!is_electron && is_future) || (is_electron && !is_future)) ? remaining : (2.0 * alpha - remaining);
    double dphi = (is_electron != is_future) ? remaining : (2.0 * alpha - remaining);
    return dphi / w;
}

ESCL_INLINE double GetBoundAngle(const double phi, const double alpha, const bool clockwise)
{
    const int multiple = (int)(phi + alpha / (PI / 2.0));
    const double bound1 = multiple * PI / 2.0 - alpha;
    const double bound2 = multiple * PI / 2.0 + alpha;

    double dangle1 = GetCrossAngle(phi, bound1, clockwise);
    double dangle2 = GetCrossAngle(phi, bound2, clockwise);

    return (dangle1 < dangle2) ? bound1 : bound2;
}

ESCL_INLINE double2 GetCyclotronOrbitCenter(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron)
{
    double2 shift = { velocity.y, -velocity.x };
    shift = shift * radius / vf;

    return is_electron ? (p - shift) : (p + shift);
}

ESCL_INLINE bool InsideImpurity(double2 pos, double2 impurity, double impurity_radius)
{
    double2 d = pos - impurity;
    return (impurity_radius * impurity_radius) > dot(d, d);
}
