#pragma once

#ifndef DEVICE_PROGRAM
	#include "src/scattering/escl/v2.h"
	#include <cmath>
#endif

#include "constants.h"
#include "details.h"

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

ESCL_INLINE int2 get_cell(const int index, const int cells_per_row)
{
	return v2i(index % cells_per_row, index / cells_per_row);
}

ESCL_INLINE int get_index(const int2 p, const int cells_per_row) {
	return p.y * cells_per_row + p.x;
}

ESCL_INLINE int2 get_cell(const double2 pos, const double2 range, const int cells_per_row)
{
	return v2i(
		(int)((pos.x - range.x) / (range.y - range.x) * (double)(cells_per_row)),
		(int)((pos.y - range.x) / (range.y - range.x) * (double)(cells_per_row))
	);
}

ESCL_INLINE double2 to_world(const int2 current_cell, const int cells_per_row, const double2 spawn_range)
{
	double factor = (spawn_range.y - spawn_range.x) / (double)cells_per_row;
	return v2(spawn_range.x + current_cell.x * factor, spawn_range.x + current_cell.y * factor);
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
	return abs(p1.x - p2.x) > 1e-6 * L || abs(p1.y - p2.y) > 1e-6 * L;
}

#define Update_Best_Intersect(intersection_point) {																			\
	double incident_angle = GetAngle(intersection_point, orbit);															\
	double dphi = GetCrossAngle(last_intersection->incident_angle, incident_angle, orbit->clockwise);									\
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

	bool horizontal_line = abs(u.x) > abs(u.y);
	bool i1_valid = horizontal_line ? PointInSegment(i1.x, p1.x, p2.x) : PointInSegment(i1.y, p1.y, p2.y);
	bool i2_valid = horizontal_line ? PointInSegment(i2.x, p1.x, p2.x) : PointInSegment(i2.y, p1.y, p2.y);

	int2 next_cell = last_intersection->entering_cell + cell_offset;

	if (i1_valid) Update_Best_Intersect(i1);
	if (i2_valid) Update_Best_Intersect(i2);

	return i1_valid || i2_valid;
}

ESCL_INLINE bool 
GetNextCell(const Orbit* orbit,
	const double phi,
	const double L,
	const int cells_per_row,
	const double2 spawn_range,
	Intersection* const last_intersection,
	Intersection* next_intersection)
{
	last_intersection->dphi = PI2;
	*next_intersection = *last_intersection;

	double2 low_left  = to_world(last_intersection->entering_cell, cells_per_row, spawn_range);
	double2 low_right = low_left + v2(L, 0);
	double2 top_right = low_left + v2(L, L);
	double2 top_left  = low_left + v2(0, L);

	UpdateFirstBoundaryIntersect(orbit, int2( 0,  1), top_left,  top_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, int2( 0, -1), low_left,  low_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, int2( 1,  0), low_right, top_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, int2(-1,  0), low_left,  top_left,  L, last_intersection, next_intersection);

	{
		double delta = 1e-7 * L;
		bool is_xlow = next_intersection->position.x - low_left.x < delta;
		bool is_xtop = -next_intersection->position.x + L + low_left.x < delta;
		bool is_ylow = next_intersection->position.y - low_left.y < delta;
		bool is_ytop = -next_intersection->position.y + L + low_left.y < delta;

		int2 offset = next_intersection->entering_cell - last_intersection->entering_cell;
		if (is_xtop && is_ylow) offset = int2(1, -1);
		if (is_xlow && is_ylow) offset = int2(-1, -1);
		if (is_xlow && is_ytop) offset = int2(-1, 1);
		if (is_xtop && is_ytop) offset = int2(1, 1);
		next_intersection->entering_cell = last_intersection->entering_cell + offset;
	}

	// Return whether we moved to a new valid cell.
	// Stop conditions:
	//	1. Progress: didn't move to a new cell
	//  2. Inside:   next cell would be outside of grid.
	//  3. Done:     next intersection would come full circle and move past the starting point.
	return (
		next_intersection->entering_cell != last_intersection->entering_cell &&
		within_bounds(next_intersection->entering_cell, cells_per_row) &&
		!AngleInRange(phi, v2(last_intersection->incident_angle, next_intersection->incident_angle), orbit->clockwise)
	);
}
