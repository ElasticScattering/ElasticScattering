#pragma once

#include "constants.h"
#include "device_macros.h"

#ifndef DEVICE_PROGRAM
    #include <cmath>
    #include <cfloat>
    #include <cstdio>
#endif

typedef struct Orbit {
	double2 center;
	double radius;
	double radius_squared;
	
	double bound_time;
	double bound_phi;

    double particle_angle;

	bool clockwise;

#ifndef DEVICE_PROGRAM
    Orbit() {};
	Orbit(double2 _center, double _radius, bool _clockwise, double _bound_time, double _bound_phi) {
		center = _center;
		radius = _radius;
		clockwise = _clockwise;
		bound_time = _bound_time;
		bound_phi = _bound_phi;

		radius_squared = radius * radius;
	}

    Orbit(double2 _center, double _radius, bool _clockwise) {
        center = _center;
        radius = _radius;
        clockwise = _clockwise;
        bound_time = 0;
        bound_phi = 0;

        radius_squared = radius * radius;
    }
#endif
} Orbit;

typedef struct Particle {
	double phi;
	double2 starting_position;
    double angular_speed;

    Orbit orbit;
} Particle;

/// 
/// Miscellaenous helper functions
/// 

ESCL_INLINE double 
smod(const double a, const double b)
{
    return a - b * floor(a / b);
}

ESCL_INLINE double 
AngleVelocity(const double pos_angle, const bool clockwise)
{
    double offset = clockwise ? HALF_PI : -HALF_PI;
    return smod(pos_angle + offset, PI2);
}


ESCL_INLINE 
double GetAngle(const double2 pos, const Orbit* orbit) {
    const double y_dist = pos.y - orbit->center.y;
    const bool left_side = pos.x < orbit->center.x;
    double angle = 0;
    if (fabs(y_dist) > (1e-10 * orbit->radius)) {
        angle = asin(y_dist / orbit->radius);
        angle = left_side ? PI - angle : angle;
    }
    else {
        angle = left_side ? PI : 0;
    }

    return AngleVelocity(angle, orbit->clockwise);
}



ESCL_INLINE double
GetCrossAngle(const double p, const double q, const bool clockwise)
{
    const double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}


ESCL_INLINE bool 
AngleInRange(const double phi, const double2 phi_range, const bool clockwise)
{
    if (fabs(phi_range.x - phi_range.y) < 1e-10) return true;

    double low  = (clockwise) ? phi_range.y : phi_range.x;
    double high = (clockwise) ? phi_range.x : phi_range.y;
    
    const double D = 1e-10;
    low  += D;
    high -= D;

    high = (high < low) ? (high + PI2) : high;
    
    return (low < phi && phi < high) || (low < (phi + PI2) && (phi + PI2) < high);
}

///
/// Impurity intersect
///
ESCL_INLINE bool 
CirclesCross(const Orbit* orbit, const double2 p2, const double r2)
{
    const double2 q = orbit->center - p2;
    const double dist_squared = dot(q, q);

    const double r_add = orbit->radius + r2;
    const double r_min = orbit->radius - r2;

    return (dist_squared < r_add * r_add) && (dist_squared > r_min * r_min);
}

ESCL_INLINE double4 
GetCrossPoints(const Orbit* orbit, const double2 p2, const double r2)
{
    const double2 q = orbit->center - p2;

    const double dist_squared = dot(q, q);
    const double dist         = sqrt(dist_squared);
    const double xs = (dist_squared + orbit->radius_squared - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(orbit->radius_squared - xs * xs);

    const double2 u = (p2 - orbit->center) / dist;

    double4 points = {
        orbit->center.x + u.x * xs +  u.y * ys,
        orbit->center.y + u.y * xs + -u.x * ys,

        orbit->center.x + u.x * xs + -u.y * ys,
        orbit->center.y + u.y * xs +  u.x * ys
    };

    return points;
}

/* Returns the first intersection with the impurity circle, or INF if it lies in an invalid range.
* 
* Orbit hits this impurity, but it can still lie in a subsection of the cell that is strictly later on 
* the orbit's trajectory than parts of other cells, and thus it should be ignored until we have found 
* no intersection in those cells.
*/
ESCL_INLINE double 
GetFirstCrossTime(const Orbit* orbit, const double2 ip, const double ir, const double w, const double2 valid_range)
{
    const double4 cross_points = GetCrossPoints(orbit, ip, ir);

    const double2 p1 = MAKE_DOUBLE2(cross_points.x, cross_points.y);
    const double2 p2 = MAKE_DOUBLE2(cross_points.z, cross_points.w);
    
    const double phi1 = GetAngle(p1, orbit);
    const double phi2 = GetAngle(p2, orbit);

    double traversal_time = INF;

    if (AngleInRange(phi1, valid_range, orbit->clockwise)) {
        const double t = GetCrossAngle(orbit->particle_angle, phi1, orbit->clockwise) / w;
        traversal_time = min(traversal_time, t);
    }

    if (AngleInRange(phi2, valid_range, orbit->clockwise)) {
        const double t = GetCrossAngle(orbit->particle_angle, phi2, orbit->clockwise) / w;
        traversal_time = min(traversal_time, t);
    }

    return traversal_time;
}

ESCL_INLINE bool
InsideImpurity(const double2 pos, const double2 impurity, const double impurity_radius)
{
    double2 d = pos - impurity;
    return (impurity_radius * impurity_radius) > dot(d, d);
}

///
/// Orbit
///

ESCL_INLINE double 
GetBoundTime(const double phi, const double alpha, const double w, const bool coherent, const bool is_electron, const bool is_future)
{
    if (coherent) return INF;

    const double remaining = smod(phi + alpha, HALF_PI);

    double dphi = (is_electron != is_future) ? remaining : (2.0 * alpha - remaining);
    return dphi / w;
}

ESCL_INLINE double 
GetBoundAngle(const double phi, const double alpha, const bool clockwise)
{
    const double v = floor((phi + alpha) / HALF_PI) * HALF_PI;
    const double bound1 = v - alpha;
    const double bound2 = v + alpha;

    double dangle1 = GetCrossAngle(phi, bound1, clockwise);
    double dangle2 = GetCrossAngle(phi, bound2, clockwise);

    return (dangle1 < dangle2) ? bound1 : bound2;
}

ESCL_INLINE double2
GetCyclotronOrbitCenter(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron)
{
    double2 shift = MAKE_DOUBLE2(velocity.y, -velocity.x );
    shift = shift * radius / vf;

    return is_electron ? (p - shift) : (p + shift);
}
