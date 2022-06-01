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

#ifndef DEVICE_PROGRAM
	#include "src/sim/es/v2.h"
	#include <cmath>
#endif

#include "constants.h"
#include "details.h"
#include "settings.h"

typedef struct Intersection {
	double2 position;
	double incident_angle;
	double dphi;
	int2 entering_cell;

#ifndef DEVICE_PROGRAM
	Intersection() {};
	Intersection(v2 pos, v2i cell, double phi)
	{
		position = pos;
		entering_cell = cell;
		incident_angle = phi;
	}
#endif
} Intersection;

ESCL_INLINE int get_index(const int2 p, const int cells_per_row) {
	return p.y * cells_per_row + p.x;
}

ESCL_INLINE int2 get_cell(const double2 pos, const double spawn_region_start, const double spawn_region_size, const int cells_per_row)
{
	return MAKE_INT2(
		(int)((pos.x - spawn_region_start) / spawn_region_size * (double)(cells_per_row)),
		(int)((pos.y - spawn_region_start) / spawn_region_size * (double)(cells_per_row))
	);
}

ESCL_INLINE double2 to_world(const int2 current_cell, const double spawn_region_start, const double cell_size)
{
	return MAKE_DOUBLE2(spawn_region_start + current_cell.x * cell_size, spawn_region_start + current_cell.y * cell_size);
}

ESCL_INLINE bool within_bounds(int2 p, const int cells_per_row) {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}

ESCL_INLINE bool PointInSegment(double point, double l0, double l1)
{
	return (point > l0 && point < l1) || (point < l0 && point > l1);
}

ESCL_INLINE bool DifferentPoint(double2 p1, double2 p2, double L)
{
	return fabs(p1.x - p2.x) > 1e-6 * L || fabs(p1.y - p2.y) > 1e-6 * L;
}


#define Update_Best_Intersect(intersection_point) {																			\
	double incident_angle = GetAngle(intersection_point, orbit);															\
	double dphi = GetCrossAngle(last_intersection->incident_angle, incident_angle, orbit->clockwise);						\
	if (dphi < closest_intersection->dphi && DifferentPoint(intersection_point, last_intersection->position, L)) {			\
		closest_intersection->position       = intersection_point;															\
		closest_intersection->entering_cell  = next_cell;																	\
		closest_intersection->incident_angle = incident_angle;																\
		closest_intersection->dphi           = dphi;																		\
	}																														\
}

ESCL_INLINE bool 
UpdateFirstBoundaryIntersect(
	const Orbit* orbit,
	const int2 cell_offset,
	const double2 p1,
	const double2 p2,
	const double L,
	const Intersection* last_intersection,
	Intersection* closest_intersection)
{
	double2 u = (p2 - p1) / L;
	double projection_distance = u.x * (orbit->center.x - p1.x) + u.y * (orbit->center.y - p1.y);
	double2 proj = p1 + (u * projection_distance);
	double proj_circle_distance_sq = pow(proj.x - orbit->center.x, 2) + pow(proj.y - orbit->center.y, 2);

	// Test if line enters circle.
	if (proj_circle_distance_sq >= orbit->radius_squared) {
		return false;
	}

	// A circle can intersect a line segment twice. Determine which 
	// intersects happened and return the earliest.
	double2 to_edge = u * sqrt(orbit->radius_squared - proj_circle_distance_sq);

	double2 i1, i2;
	i1 = proj + to_edge;
	i2 = proj - to_edge;

	bool horizontal_line = fabs(u.x) > fabs(u.y);
	bool i1_valid = horizontal_line ? PointInSegment(i1.x, p1.x, p2.x) : PointInSegment(i1.y, p1.y, p2.y);
	bool i2_valid = horizontal_line ? PointInSegment(i2.x, p1.x, p2.x) : PointInSegment(i2.y, p1.y, p2.y);

	int2 next_cell = last_intersection->entering_cell + cell_offset;

	if (i1_valid) Update_Best_Intersect(i1);
	
	if (i2_valid) Update_Best_Intersect(i2);
	return i1_valid || i2_valid;
}


ESCL_INLINE bool 
GetNextCell(
	const Orbit* orbit,
	const double phi,
	IMPURITY_SETTINGS,
	Intersection* const last_intersection,
	Intersection* next_intersection)
{
	last_intersection->dphi = PI2;
	*next_intersection = *last_intersection;

	const double L = is->cell_size;
	double2 low_left  = to_world(last_intersection->entering_cell, is->spawn_region_start, L);
	double2 low_right = low_left + MAKE_DOUBLE2(L, 0);
	double2 top_right = low_left + MAKE_DOUBLE2(L, L);
	double2 top_left  = low_left + MAKE_DOUBLE2(0, L);

	UpdateFirstBoundaryIntersect(orbit, MAKE_INT2( 0,  1), top_left,  top_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, MAKE_INT2( 0, -1), low_left,  low_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, MAKE_INT2( 1,  0), low_right, top_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, MAKE_INT2(-1,  0), low_left,  top_left,  L, last_intersection, next_intersection);

	{
		double delta = 1e-7 * L;
		bool is_xlow = ( next_intersection->position.x     - low_left.x) < delta;
		bool is_xtop = (-next_intersection->position.x + L + low_left.x) < delta;
		bool is_ylow = ( next_intersection->position.y     - low_left.y) < delta;
		bool is_ytop = (-next_intersection->position.y + L + low_left.y) < delta;

		int2 offset = next_intersection->entering_cell - last_intersection->entering_cell;
		if (is_xtop && is_ylow) offset = MAKE_INT2( 1, -1);
		if (is_xlow && is_ylow) offset = MAKE_INT2(-1, -1);
		if (is_xlow && is_ytop) offset = MAKE_INT2(-1,  1);
		if (is_xtop && is_ytop) offset = MAKE_INT2( 1,  1);
		next_intersection->entering_cell = last_intersection->entering_cell + offset;
	}

	// Return whether we moved to a new valid cell.
	// Stop conditions:
	//	1. Progress: next cell is the same as current cell.
	//  2. Inside:   next cell would be outside of grid.
	//  3. Done:     next intersection would come full circle and move past the starting point.
	#define SAME_CELL(c1, c2) c1.x == c2.x && c1.y == c2.y
	bool same_cell = SAME_CELL(next_intersection->entering_cell, last_intersection->entering_cell);

	return
		!same_cell &&
		within_bounds(next_intersection->entering_cell, is->cells_per_row) &&
		!AngleInRange(phi, MAKE_DOUBLE2(last_intersection->incident_angle, next_intersection->incident_angle), orbit->clockwise);
}
