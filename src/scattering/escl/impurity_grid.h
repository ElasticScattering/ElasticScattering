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
	int entering_cell;
} Intersection;

struct CellRange {
	int start, end;
};

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
	return (point > l0 && point < l1) || (point < l0&& point > l1);
}

inline bool GetFirstBoundaryIntersect(const double2 p1, const double2 p2, const Orbit* orbit, const double L, const double start_phi, Intersection* intersection) {
	double2 u = (p2 - p1) / L;
	double projection_distance = u.x * (orbit->center.x - p1.x) + u.y * (orbit->center.y - p1.y);
	double2 proj = p1 + (u * projection_distance);
	double proj_circle_distance_sq = pow(proj.x - orbit->center.x, 2) + pow(proj.y - orbit->center.y, 2);

	// Test if line enters circle.
	if (proj_circle_distance_sq >= orbit->radius_squared) {
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

inline bool GetNextCell(const Orbit* orbit,
	const double2 current_cell_lowleft,
	const Intersection last_intersection,
	const double L,
	const int cells_per_row,
	Intersection* next_intersection)
{
	double2 low_left = current_cell_lowleft;
	double2 low_right = low_left + double2(L, 0);
	double2 top_right = low_left + double2(L, L);
	double2 top_left = low_left + double2(0, L);

	Intersection left, right, up, down;

	double last_phi = last_intersection.dphi;
	bool up_hit = GetFirstBoundaryIntersect(top_left, top_right, orbit, L, last_phi, &up);
	bool down_hit = GetFirstBoundaryIntersect(low_left, low_right, orbit, L, last_phi, &down);
	bool right_hit = GetFirstBoundaryIntersect(low_right, top_right, orbit, L, last_phi, &right);
	bool left_hit = GetFirstBoundaryIntersect(low_left, top_left, orbit, L, last_phi, &left);

	if (!up_hit && !down_hit && !right_hit && !left_hit)
		return false;

	// Select one valid intersect.
	Intersection closest_intersection = last_intersection;

	if (up_hit) {
		// Test if this intersection is closer than what we have
		double dphi = GetCrossAngle(last_intersection.incident_angle, up.incident_angle, orbit->clockwise);
		int next_cell_candidate = last_intersection.entering_cell - cells_per_row;
		if (dphi < closest_intersection.dphi && next_cell_candidate >= 0) {
			closest_intersection.dphi = dphi;
			closest_intersection.entering_cell = last_intersection.entering_cell - cells_per_row;
		}
	}

	if (down_hit) {
		// Test if this intersection is closer than what we have
		double dphi = GetCrossAngle(last_intersection.incident_angle, up.incident_angle, orbit->clockwise);
		int next_cell_candidate = last_intersection.entering_cell + cells_per_row; // @Todo, dit controleren..
		if (dphi < closest_intersection.dphi && next_cell_candidate >= 0) {
			closest_intersection.dphi = dphi;
			closest_intersection.entering_cell = last_intersection.entering_cell - cells_per_row; // ....
		}
	}
	//etc ...

	if (closest_intersection.entering_cell == last_intersection.entering_cell) {
		return false;
	}

	// The intersection can not put us out of bounds.
	//@Todo: doe dit met cell_x, cell_y, tijdens de 4-test, of achteraf?

	return true;
}

inline v2i to_grid(const double x, const double y, const double2 range, const int cells_per_row)
{
	return {
		(int)((x - range.x) / (range.y - range.x) * (cells_per_row)),
		(int)((y - range.x) / (range.y - range.x) * (cells_per_row))
	};
}

inline v2 to_world(const int cell_index, const int cells_per_row, const double2 spawn_range)
{

	double x = cell_index % cells_per_row;
	double y = floor(cell_index / cells_per_row);
	double2 low_left = (double2)(x, y) * (spawn_range.y - spawn_range.x) + spawn_range.x;

	return low_left;
}

inline bool within_bounds(v2i p, const int cells_per_row) {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}

inline int to_index(const v2i p, const int cells_per_row) {
	return p.y * cells_per_row + p.x;
}

////////////////////////////

inline int get_cell_index(const v2 pos, const v2 range, const int cells_per_row)
{
	return to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
}

inline CellRange get_cell_id(const v2 pos, const v2 range, const int cells_per_row)
{
	CellRange cell_range;

	cell_range.start = to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
	cell_range.end = cell_range.start + 1;

	return cell_range;
}

////////////////////////

inline double remap(double x, double s0, double s1, double t0, double t1)
{
	return t0 + (x - s0) / (s1 - s0) * (t1 - t0);
}
