#pragma once

#include "v2.h"

struct Orbit {
	double2 position;
	double radius;
	double radius_squared;

	Orbit(double2 pos, double r) {
		position = pos;
		radius = r;
		radius_squared = r * r;
	}
};
