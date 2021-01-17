#pragma once

#include "src/sim/es/shared_macros.h"

#ifdef DEVICE_PROGRAM
	#define IMPURITY_SETTINGS __constant ImpuritySettings *is
	#define PARTICLE_SETTINGS __constant ParticleSettings* ps

	#define BUFFER_ARGS __constant double2 *impurities, __constant int *cell_indices, __global ParticleMetrics *metrics

	#define MAKE_DOUBLE2(px, py) (double2)(px, py)
	#define MAKE_INT2(px, py) (int2)(px, py)

	#define METRIC_INC_COND(cond, metric_name) ((cond)) atomic_inc(&metrics->metric_name)
	#define METRIC_INC(metric_name)						atomic_inc(&metrics->metric_name);

	#define METRIC_ADD(metric_name, m_value)			atomic_add(&metrics->metric_name, (m_value))

#else
	#include "v2.h"

	#define IMPURITY_SETTINGS const ImpuritySettings* const is
	#define PARTICLE_SETTINGS const ParticleSettings* const ps
	#define BUFFER_ARGS const std::vector<v2> &impurities, const std::vector<int> &cell_indices, ParticleMetrics &metrics

	#define MAKE_DOUBLE2(px, py) v2(px, py)
	#define MAKE_INT2(px, py) v2i(px, py)

	#define min(a, b) ((a) < (b)) ? (a) : (b)
	#define max(a, b) ((a) > (b)) ? (a) : (b)

	#define METRIC_INC(metric_name) metrics.metric_name += 1
	#define METRIC_INC_COND(cond, metric_name) if ((cond)) metrics.metric_name += 1
	#define METRIC_ADD(metric_name, m_value) metrics.metric_name += (m_value)
#endif
