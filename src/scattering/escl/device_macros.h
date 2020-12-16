#pragma once

#ifdef DEVICE_PROGRAM
	#define BUFFER_ARGS __constant ScatteringParameters* sp, __global double2* impurities, __global int* cell_indices
	

	#define MAKE_ORBIT(center, r, clockwise, bound_time, bound_phi) (Orbit){ .center = center, .radius = r, .clockwise = clockwise, .bound_time = bound_time, .bound_phi = bound_phi, .radius_squared = (r*r)}
	#define MAKE_POS(i, j) MAKE_DOUBLE(i, j) * (sp.region_size / (double)(sp.dim - 2))	

	#define MAKE_DOUBLE2(px, py) (double2)(px, py)
	#define MAKE_INT2(px, py) (int2)(px, py)

	#define ESCL_INLINE
#else
	#include "v2.h"

	#define BUFFER_ARGS const ScatteringParameters *sp, const std::vector<v2> &impurities, const std::vector<int> &cell_indices, Metrics *metrics

	#define MAKE_ORBIT(center, radius, clockwise, bound_time, bound_phi) Orbit(center, radius, clockwise, bound_time, bound_phi)
	#define MAKE_DOUBLE2(px, py) v2(px, py)
	#define MAKE_INT2(px, py) v2i(px, py)

	#define ESCL_INLINE inline

	#define min(a, b) ((a) < (b)) ? (a) : (b)
	#define max(a, b) ((a) > (b)) ? (a) : (b)
#endif
