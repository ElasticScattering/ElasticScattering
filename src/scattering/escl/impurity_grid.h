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
} Intersection;

struct CellRange {
	int start, end;
};


inline int2 get_cell(const int index, const int cells_per_row);
inline int get_index(const int2 p, const int cells_per_row);
inline bool within_bounds(int2 p, const int cells_per_row);
inline double2 to_world(const int2 current_cell, const int cells_per_row, const double2 spawn_range);

inline double GetAngle(double2 pos, const Orbit* orbit) {
	if (abs(pos.y - orbit->center.y) > EPSILON) {
		double sign = (pos.x < orbit->center.x) ? -1.0 : 1.0;
		double angle = sign * asin((pos.y - orbit->center.y) / orbit->radius);
		return fmod(angle, PI2);
	}
	else {
		return (pos.x > orbit->center.x) ? 0 : PI;
	}
}

inline bool PointInSegment(double point, double l0, double l1)
{
	return (point > l0 && point < l1) || (point < l0 && point > l1);
}

inline bool GetFirstBoundaryIntersect(const double2 p1, const double2 p2, const double L, const Orbit* orbit, const double start_phi, Intersection* intersection) {
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

	double2 i1 = proj + to_edge;
	double2 i2 = proj - to_edge;

	// Determine if the intersection happen on the line segment.
	bool horizontal_line = abs(u.x) > abs(u.y);

	Intersection in1;
	in1.position = i1;
	in1.incident_angle = GetAngle(i1, orbit);
	in1.dphi = GetCrossAngle(start_phi, in1.incident_angle, orbit->clockwise);

	Intersection in2;
	in2.position = i2;
	in2.incident_angle = GetAngle(i2, orbit);
	in2.dphi = GetCrossAngle(start_phi, in2.incident_angle, orbit->clockwise);

	bool i1_valid = horizontal_line ? PointInSegment(i1.x, p1.x, p2.x) : PointInSegment(i1.y, p1.y, p2.y);
	bool i2_valid = horizontal_line ? PointInSegment(i2.x, p1.x, p2.x) : PointInSegment(i2.y, p1.y, p2.y);

	//@Refactor, kan dit simpeler omdat clockwise al in GetCrossAngle zit?
	if (i1_valid && i2_valid)
	{
		bool phi1_lower = in1.dphi < in2.dphi;
		if (orbit->clockwise) intersection = (phi1_lower) ? &in1 : &in2;
		else             intersection = (phi1_lower) ? &in2 : &in1;
	}
	else if (i1_valid) intersection = &in1;
	else               intersection = &in2;

	return true;
}

inline void WriteBestIntersect(const Intersection last_intersection, int2 offset, double incident_angle, bool clockwise, int cells_per_row, Intersection *closest_intersection)
{
	int2 next_cell_candidate = last_intersection.entering_cell + offset;

	double dphi = GetCrossAngle(last_intersection.incident_angle, incident_angle, clockwise);
	if (dphi < closest_intersection->dphi && within_bounds(next_cell_candidate, cells_per_row)) {
		closest_intersection->dphi = dphi;
		closest_intersection->entering_cell = next_cell_candidate;
	}
}

inline bool GetNextCell(const Orbit* orbit,
	const Intersection last_intersection,
	const double L,
	const int cells_per_row,
	const double2 spawn_range,
	Intersection* next_intersection)
{
	double2 low_left  = to_world(last_intersection.entering_cell, cells_per_row, spawn_range);
	double2 low_right = low_left + double2(L, 0);
	double2 top_right = low_left + double2(L, L);
	double2 top_left  = low_left + double2(0, L);

	Intersection left, right, up, down;
	bool up_hit    = GetFirstBoundaryIntersect(top_left , top_right, L, orbit, last_intersection.dphi, &up);
	bool down_hit  = GetFirstBoundaryIntersect(low_left , low_right, L, orbit, last_intersection.dphi, &down);
	bool right_hit = GetFirstBoundaryIntersect(low_right, top_right, L, orbit, last_intersection.dphi, &right);
	bool left_hit  = GetFirstBoundaryIntersect(low_left , top_left , L, orbit, last_intersection.dphi, &left);

	Intersection closest = last_intersection;
	if (up_hit)    WriteBestIntersect(last_intersection, int2(0, 1) , up.incident_angle,    orbit->clockwise, cells_per_row, &closest);
	if (down_hit)  WriteBestIntersect(last_intersection, int2(0, -1), down.incident_angle,  orbit->clockwise, cells_per_row, &closest);
	if (right_hit) WriteBestIntersect(last_intersection, int2(1, 0) , right.incident_angle, orbit->clockwise, cells_per_row, &closest);
	if (left_hit)  WriteBestIntersect(last_intersection, int2(-1, 0), left.incident_angle,  orbit->clockwise, cells_per_row, &closest);

	// Return whether we moved to a new cell.
	*next_intersection = closest;
	return (next_intersection->entering_cell != last_intersection.entering_cell);
}

inline int2 get_cell(const double x, const double y, const double2 range, const int cells_per_row)
{
	return {
		(int)((x - range.x) / (range.y - range.x) * (double)(cells_per_row)),
		(int)((y - range.x) / (range.y - range.x) * (double)(cells_per_row))
	};
}

inline int2 get_cell(const int index, const int cells_per_row)
{
	return (int2)(index % cells_per_row, index / cells_per_row);
}

inline int get_index(const int2 p, const int cells_per_row) {
	return p.y * cells_per_row + p.x;
}

inline int get_cell_index(const double2 pos, const double2 range, const int cells_per_row)
{
	return get_index(get_cell(pos.x, pos.y, range, cells_per_row), cells_per_row);
}

inline double2 to_world(const int2 current_cell, const int cells_per_row, const double2 spawn_range)
{
	return (double2)(current_cell.x, current_cell.y) * (spawn_range.y - spawn_range.x) + (double2)(spawn_range.x, spawn_range.x);
}

inline bool within_bounds(int2 p, const int cells_per_row) {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}

/*
inline double remap(double x, double s0, double s1, double t0, double t1)
{
	return t0 + (x - s0) / (s1 - s0) * (t1 - t0);
}
*/
