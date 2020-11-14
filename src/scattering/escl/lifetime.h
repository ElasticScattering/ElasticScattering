#pragma once

#include "ScatteringParameters.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "impurity_grid.h"

#include <vector>
#include "windows.h"


double lifetime(const int quadrant, const int step, const double2 pos, BUFFER_ARGS);
double SingleLifetime(const Particle* p, const double2 impurity, const double impurity_radius, const double angular_speed);
double TraceOrbit(Particle* p, BUFFER_ARGS);
