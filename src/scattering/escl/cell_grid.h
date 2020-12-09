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

struct CellRange {
	int start, end;
};


ESCL_INLINE int2 get_cell(const int index, const int cells_per_row);
ESCL_INLINE int get_index(const int2 p, const int cells_per_row);
ESCL_INLINE bool within_bounds(int2 p, const int cells_per_row);
ESCL_INLINE double2 to_world(const int2 current_cell, const int cells_per_row, const double2 spawn_range);

ESCL_INLINE bool PointInSegment(double point, double l0, double l1)
{
	return (point > l0 && point < l1) || (point < l0 && point > l1);
}

ESCL_INLINE bool DifferentPoint(double2 p1, double2 p2, double L)
{
	return abs(p1.x - p2.x) > 1e-6 * L || abs(p1.y - p2.y) > 1e-6 * L;
}

#if 0
ESCL_INLINE void UpdateBestIntersectOld(Intersection candidate, int2 offset, const Intersection *last_intersection, bool clockwise, double L, Intersection* closest_intersection)
{
	candidate.entering_cell = last_intersection->entering_cell + offset;
	candidate.dphi = GetCrossAngle(last_intersection->incident_angle, candidate.incident_angle, clockwise);

	//printf("\tComparing Intersection: (%e, %e, %e)\n", candidate.position.x, candidate.position.y, candidate.dphi);
	//printf("\tComparing to: (%e, %e, %e)\n", closest_intersection->position.x, closest_intersection->position.y, closest_intersection->dphi);

	if (candidate.dphi < closest_intersection->dphi && 
		DifferentPoint(candidate.position, last_intersection->position, L)) {
		*closest_intersection = candidate;
	}
}

ESCL_INLINE bool GetFirstBoundaryIntersectOld(const Orbit* orbit, const double2 p1, const double2 p2, const double L, const double start_phi, Intersection* intersection) {
	double2 u = (p2 - p1) / L;
	double projection_distance = u.x * (orbit->center.x - p1.x) + u.y * (orbit->center.y - p1.y);
	double2 proj = p1 + (u * projection_distance);
	double proj_circle_distance_sq = pow(proj.x - orbit->center.x, 2) + pow(proj.y - orbit->center.y, 2);

	// Test if line enters circle.
	if (proj_circle_distance_sq >= orbit->radius_squared) {
		intersection->dphi = INF; // Om te checken.
		return false;
	}

	// Calculate segment to the edge of the circle.
	double2 to_edge = u * sqrt(orbit->radius_squared - proj_circle_distance_sq);

	// A circle can intersect a line segment twice. Determine which 
	// intersects happened and return the earliest.
	Intersection i1, i2;
	i1.position = proj + to_edge;
	i2.position = proj - to_edge;

	bool horizontal_line = abs(u.x) > abs(u.y);
	bool i1_valid = horizontal_line ? PointInSegment(i1.position.x, p1.x, p2.x) : PointInSegment(i1.position.y, p1.y, p2.y);
	bool i2_valid = horizontal_line ? PointInSegment(i2.position.x, p1.x, p2.x) : PointInSegment(i2.position.y, p1.y, p2.y);

	if (!i1_valid && !i2_valid) return false;

	i1.incident_angle = GetAngle(i1.position, orbit);
	i1.dphi = GetCrossAngle(start_phi, i1.incident_angle, orbit->clockwise);

	i2.incident_angle = GetAngle(i2.position, orbit);
	i2.dphi = GetCrossAngle(start_phi, i2.incident_angle, orbit->clockwise);

	//@Refactor, kan dit simpeler omdat clockwise al in GetCrossAngle zit?
	if (i1_valid && i2_valid)
	{
		bool phi1_lower = i1.dphi < i2.dphi;
		if (orbit->clockwise) *intersection = (phi1_lower) ? i1 : i2;
		else				  *intersection = (phi1_lower) ? i2 : i1;
	}
	else if (i1_valid) *intersection = i1;
	else               *intersection = i2;

	return true;
}

ESCL_INLINE bool GetNextCellOld(const Orbit* orbit,
	Intersection * const last_intersection,
	const double L,
	const int cells_per_row,
	const double2 spawn_range,
	Intersection* next_intersection)
{
	double2 low_left  = to_world(last_intersection->entering_cell, cells_per_row, spawn_range);
	double2 low_right = low_left + v2(L, 0);
	double2 top_right = low_left + v2(L, L);
	double2 top_left  = low_left + v2(0, L);

	printf("Cell: (%i, %i) -> Low left: (%e, %e)\n", last_intersection->entering_cell.x, last_intersection->entering_cell.y, low_left.x, low_left.y);

	last_intersection->dphi = PI2;

	Intersection i_left, i_right, i_up, i_down;
	bool hit_up    = GetFirstBoundaryIntersectOld(orbit, top_left,  top_right, L, last_intersection->dphi, &i_up);
	bool hit_down  = GetFirstBoundaryIntersectOld(orbit, low_left,  low_right, L, last_intersection->dphi, &i_down);
	bool hit_right = GetFirstBoundaryIntersectOld(orbit, low_right, top_right, L, last_intersection->dphi, &i_right);
	bool hit_left  = GetFirstBoundaryIntersectOld(orbit, low_left,  top_left,  L, last_intersection->dphi, &i_left);

	*next_intersection = *last_intersection;
	if (hit_up)    UpdateBestIntersectOld(i_up,    int2(0,  1), last_intersection, orbit->clockwise, L, next_intersection);
	if (hit_down)  UpdateBestIntersectOld(i_down,  int2(0, -1), last_intersection, orbit->clockwise, L, next_intersection);
	if (hit_right) UpdateBestIntersectOld(i_right, int2(1, 0),  last_intersection, orbit->clockwise, L, next_intersection);
	if (hit_left)  UpdateBestIntersectOld(i_left,  int2(-1, 0), last_intersection, orbit->clockwise, L, next_intersection);
	
	//printf("Closest Intersection: (%e, %e)\n", next_intersection->position.x, next_intersection->position.y);
	//printf("\t In bounds: " + ((within_bounds(candidate.entering_cell, cells_per_row)) ? "True" : "False") + "\n");

	// Return whether we moved to a new valid cell.
	// If we didn't move or to a cell that is out of bounds the search stops.
	return (next_intersection->entering_cell != last_intersection->entering_cell) && within_bounds(next_intersection->entering_cell, cells_per_row);
}
#endif

ESCL_INLINE void UpdateBestIntersect(const int2 next_cell, const double2 intersection_point, const Intersection* last_intersection, const Orbit* orbit, double L, Intersection* closest_intersection)
{
	Intersection candidate;
	candidate.position = intersection_point;
	candidate.entering_cell = next_cell;
	candidate.incident_angle = GetAngle(candidate.position, orbit);
	candidate.dphi = GetCrossAngle(last_intersection->dphi, candidate.incident_angle, orbit->clockwise);

	if (candidate.dphi < closest_intersection->dphi &&
		DifferentPoint(candidate.position, last_intersection->position, L)) {
		
		//printf("\t\tPhi Candidate: %f :: (%e, %e) \n", candidate.incident_angle, candidate.position.x, candidate.position.y);

		*closest_intersection = candidate;
	}
}

ESCL_INLINE bool UpdateFirstBoundaryIntersect(
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
		closest_intersection->dphi = INF; // Om te checken.
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

	if (i1_valid) UpdateBestIntersect(next_cell, i1, last_intersection, orbit, L, closest_intersection);
	if (i2_valid) UpdateBestIntersect(next_cell, i2, last_intersection, orbit, L, closest_intersection);

	return i1_valid || i2_valid;
}


ESCL_INLINE bool GetNextCell(const Orbit* orbit,
	const double phi,
	Intersection* const last_intersection,
	const double L,
	const int cells_per_row,
	const double2 spawn_range,
	Intersection* next_intersection)
{
	last_intersection->dphi = PI2;
	*next_intersection = *last_intersection;

	double2 low_left  = to_world(last_intersection->entering_cell, cells_per_row, spawn_range);
	double2 low_right = low_left + v2(L, 0);
	double2 top_right = low_left + v2(L, L);
	double2 top_left  = low_left + v2(0, L);

	bool any_intersections = false;
	UpdateFirstBoundaryIntersect(orbit, int2( 0,  1), top_left,  top_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, int2( 0, -1), low_left,  low_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, int2( 1,  0), low_right, top_right, L, last_intersection, next_intersection);
	UpdateFirstBoundaryIntersect(orbit, int2(-1,  0), low_left,  top_left,  L, last_intersection, next_intersection);

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

ESCL_INLINE int2 get_cell(const double x, const double y, const double2 range, const int cells_per_row)
{
	return v2i(
		(int)((x - range.x) / (range.y - range.x) * (double)(cells_per_row)),
		(int)((y - range.x) / (range.y - range.x) * (double)(cells_per_row))
	);
}

ESCL_INLINE int2 get_cell(const int index, const int cells_per_row)
{
	return v2i(index % cells_per_row, index / cells_per_row);
}

ESCL_INLINE int get_index(const int2 p, const int cells_per_row) {
	return p.y * cells_per_row + p.x;
}

ESCL_INLINE int get_cell_index(const double2 pos, const double2 range, const int cells_per_row)
{
	return get_index(get_cell(pos.x, pos.y, range, cells_per_row), cells_per_row);
}

ESCL_INLINE double2 to_world(const int2 current_cell, const int cells_per_row, const double2 spawn_range)
{
	double factor = (spawn_range.y - spawn_range.x) / (double)cells_per_row;
	return v2(spawn_range.x + current_cell.x * factor, spawn_range.x + current_cell.y * factor);
}

ESCL_INLINE bool within_bounds(int2 p, const int cells_per_row) {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}
