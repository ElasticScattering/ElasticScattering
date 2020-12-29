#pragma once

#ifdef DEVICE_PROGRAM
	#define ESCL_INLINE

#define IMPURITY_SETTINGS __constant ImpuritySettings *settings
	#define PARTICLE_SETTINGS __constant ParticleSettings* ps

	#define BUFFER_ARGS __global double2 *impurities, __global int *cell_indices, __global Metrics *metrics
	
	#define MAKE_POS(i, j) MAKE_DOUBLE(i, j) * (sp.region_size / (double)(sp.dim - 2))	

	#define MAKE_DOUBLE2(px, py) (double2)(px, py)
	#define MAKE_INT2(px, py) (int2)(px, py)

	#define METRIC_INC_COND(cond, metric_name) if ((cond)) atomic_inc(&metric_name)
	#define METRIC_INC(metric_name) atomic_inc(&metric_name)
	#define METRIC_ADD(metric_name, m_value) atomic_add(&metric_name, (m_value))

#else
	#include "v2.h"

	#define ESCL_INLINE inline

	#define IMPURITY_SETTINGS const ImpuritySettings * const settings
	#define PARTICLE_SETTINGS const ParticleSettings* const ps
	#define BUFFER_ARGS const std::vector<v2> &impurities, const std::vector<int> &cell_indices, Metrics *metrics

	#define MAKE_DOUBLE2(px, py) v2(px, py)
	#define MAKE_INT2(px, py) v2i(px, py)

	#define min(a, b) ((a) < (b)) ? (a) : (b)
	#define max(a, b) ((a) > (b)) ? (a) : (b)

	#define METRIC_INC(metric_name) metric_name += 1
	#define METRIC_INC_COND(cond, metric_name) if ((cond)) metric_name += 1
	#define METRIC_ADD(metric_name, m_value) metric_name += (m_value)
#endif
