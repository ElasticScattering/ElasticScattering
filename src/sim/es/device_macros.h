/* -------------------------------------------------------------------------
	This code is part of ElasticScattering.

	Copyright(C) 2022 Elastic Scattering developers

	This program is free software : you can redistribute it and /or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

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
