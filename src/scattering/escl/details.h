#pragma once
#include "constants.h"
#include "device_macros.h"

typedef struct Orbit {
	double2 center;
	double radius;
	double radius_squared;
	
	bool clockwise;

	Orbit(double2 _center, double _radius, bool _clockwise) {
		center = _center;
		radius = _radius;
		clockwise = _clockwise;

		radius_squared = radius * radius;
	}
} Orbit;

double smod(const double a, const double b);

double GetBoundTime(const double phi, const double alpha, const double w, const bool is_incoherent, const bool is_diag_region, const bool is_electron, const bool is_future);
double2 GetCyclotronOrbitCenter(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron);

bool CirclesCross(const Orbit c1, const double2 p2, const double r2);
double4 GetCrossPoints(const Orbit c1, const double2 p2, const double r2);
double GetFirstCrossTime(const double2 pos, const Orbit c1, const double2 ip, const double ir, const double w);

double GetPhi(const double2 pos, const Orbit c1);
double GetCrossAngle(const double p, const double q, const bool clockwise);
