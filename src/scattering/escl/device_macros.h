#pragma once

#ifdef DEVICE_PROGRAM
	#define BUFFER_ARGS __constant ScatteringParameters* sp, __global read_only double2* impurities, __global read_only int* cell_indices
	#define BUFFER_ARGS3D __constant ScatteringParameters* sp, __global read_only double2* impurities, __local double* phi_lt
#else
	#include "v2.h"

	#define BUFFER_ARGS const ScatteringParameters *sp, const std::vector<v2> &impurities, const std::vector<int> &cell_indices
	#define BUFFER_ARGS3D const ScatteringParameters *sp, const std::vector<v2> &impurities, const std::vector<double> phi_lt

	#define min(a, b) ((a) < (b)) ? (a) : (b)
#endif
