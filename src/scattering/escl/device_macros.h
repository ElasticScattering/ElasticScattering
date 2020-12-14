#pragma once

#ifdef DEVICE_PROGRAM
	#define BUFFER_ARGS __constant ScatteringParameters* sp, __global read_only double2* impurities, __global read_only int* cell_indices
	
	#define ESCL_INLINE
#else
	#include "v2.h"

	#define BUFFER_ARGS const ScatteringParameters *sp, const std::vector<v2> &impurities, const std::vector<int> &cell_indices, Metrics *metrics

	#define ESCL_INLINE inline

	#define min(a, b) ((a) < (b)) ? (a) : (b)
	#define max(a, b) ((a) > (b)) ? (a) : (b)

#endif
