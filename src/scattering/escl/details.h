#pragma once
#include "constants.h"
#include "device_macros.h"

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
	Orbit orbit;
	int cell_index;
	double phi;
	double2 starting_position;
} Particle;

double smod(const double a, const double b);

double GetBoundTime(const double phi, const double alpha, const double w, const bool is_incoherent, const bool is_diag_region, const bool is_electron, const bool is_future);
double GetBoundAngle(const double phi, const double alpha, const bool clockwise);
double2 GetCyclotronOrbitCenter(const double2 p, const double2 velocity, const double radius, const double vf, const bool is_electron);

bool CirclesCross(const Orbit c1, const double2 p2, const double r2);
double4 GetCrossPoints(const Orbit c1, const double2 p2, const double r2);
double GetFirstCrossTime(const double2 pos, const Orbit c1, const double2 ip, const double ir, const double w);

double GetPhi(const double2 pos, const Orbit c1);
double GetCrossAngle(const double p, const double q, const bool clockwise);
