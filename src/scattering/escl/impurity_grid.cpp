#pragma once

#include "impurity_grid.h"

#ifndef DEVICE_PROGRAM
	#include <cmath>
#endif
#include "constants.h"

double GetAngle(double2 pos, Orbit o) {
	if (abs(pos.y - o.center.y) > EPSILON) {
		double sign = (pos.x < o.center.x) ? -1.0 : 1.0;
		double angle = sign * asin((pos.y - o.center.y) / o.radius);
		return fmod(angle, PI2);
	}
	else {
		return (pos.x > o.center.x) ? 0 : PI;
	}
}

bool PointInSegment(double point, double l0, double l1)
{
	return ((point > l0) && (point < l1)) || ((point < l0) && (point > l1));
}

bool GetBoundaryIntersects(const double2 p1, const double2 p2, const Orbit o, const double L, SideIntersection *intersection) {
	double2 u = (p2 - p1) / L;
	double projection_distance = u.x * (o.center.x - p1.x) + u.y * (o.center.y - p1.y);
	double2 proj = p1 + (u * projection_distance);
	double proj_circle_distance_sq = pow(proj.x - o.center.x, 2) + pow(proj.y - o.center.y, 2);
	
	// Test if line enters circle.
	if (proj_circle_distance_sq >= o.radius_squared) {
		return false; 
	}
	
	// Calculate segment to the edge of the circle.
	double2 to_edge = u * sqrt(o.radius_squared - proj_circle_distance_sq);

	double2 i1 = proj + to_edge;
	double2 i2 = proj - to_edge;

	// Determine if the intersection happen on the line segment.
	bool horizontal_line = abs(u.x) > abs(u.y);

	bool i1_valid = horizontal_line ? PointInSegment(i1.x, p1.x, p2.x) : PointInSegment(i1.y, p1.y, p2.y);
	bool i2_valid = horizontal_line ? PointInSegment(i2.x, p1.x, p2.x) : PointInSegment(i2.y, p1.y, p2.y);

	// @Todo, decide which to return?
	if (i1_valid) {
		Intersection i;
		i.position = i1;
		i.incident_angle = GetAngle(i.position, o);
		i.dphi = GetCrossAngle(, i.incident_angle, o.clockwise);
		intersection->i1 = i;
	}
	if (i2_valid) {
		Intersection i;
		i.position = i2;
		i.incident_angle = GetAngle(i.position, o);
		i.dphi = GetCrossAngle(, i.incident_angle, o.clockwise);
		intersection->i2 = i;
	}

	return true;
}

void FindExitIntersect(Orbit o, double phi, Intersection start_intersect)
{
	// IN: Intersects ...

	double dphi_closest = GetCrossAngle(start_intersect.incident_angle, phi, o.clockwise);

// for intersects
	Intersection i; // from for loop
	double dphi = GetCrossAngle(start_intersect.incident_angle, i.incident_angle, o.clockwise);
	if (dphi < dphi_closest) {
		// save this intersect?
	}
}

bool GetNextCell(
	const int current_cell,
	const Orbit orbit, 
	const double2 last_intersection,
	const int cells_per_row,
	const double L,
	const double2 spawn_range,
	int *next_cell,
	double2* intersection_point)
{
	double2 low_left = to_world(current_cell, cells_per_row, spawn_range);
	double2 top_right = low_left + double2(L, L);
	double2 top_left = low_left + double2(0, L);
	double2 low_right = low_left + double2(L, 0);

	SideIntersection up, down, left, right;
	bool up_hit    = GetBoundaryIntersects(low_left, top_left, orbit, L, &up);
	bool down_hit  = GetBoundaryIntersects(low_right, top_right, orbit, L, &down);
	bool left_hit  = GetBoundaryIntersects(top_left, top_right, orbit, L, &left);
	bool right_hit = GetBoundaryIntersects(low_left, low_right, orbit, L, &right);

	if (!right_hit && !up_hit && !left_hit && !down_hit)
		return false;

	// Select [one] from all valid intersects.
	// Some succeeded, with 1 or 2 intersects....

	// One attempt would be 

	if (up_hit) {
		// One or two?
	}

	return true;
}

////////////////////////////

int get_cell_index(const v2 pos, const v2 range, const int cells_per_row)
{
	return to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
}

CellRange get_cell_id(const v2 pos, const v2 range, const int cells_per_row)
{
	CellRange cell_range;

	cell_range.start = to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
	cell_range.end = cell_range.start + 1;

	return cell_range;
}

////////////////////////


int to_grid(const double x, const double2 range, const int cells_per_row)
{
	return (int)((x - range.x) / (range.y - range.x) * (cells_per_row));
}

v2i to_grid(const double x, const double y, const double2 range, const int cells_per_row)
{
	return {
		(int)((x - range.x) / (range.y - range.x) * (cells_per_row)),
		(int)((y - range.x) / (range.y - range.x) * (cells_per_row))
	};
}

v2 to_world(const int cell_index, const int cells_per_row, const double2 spawn_range)
{

	double x = cell_index % cells_per_row;
	double y = floor(cell_index / cells_per_row);
	double2 low_left = (double2)(x, y) * (spawn_range.y - spawn_range.x) + spawn_range.x;

	return low_left;
}

bool within_bounds(v2i p, const int cells_per_row) {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}

int to_index(const int x, const int y, const int cells_per_row) {
	return y * cells_per_row + x;
}

int to_index(const v2i p, const int cells_per_row) {
	return to_index(p.x, p.y, cells_per_row);
}

double remap(double x, double s0, double s1, double t0, double t1)
{
	return t0 + (x - s0) / (s1 - s0) * (t1 - t0);
}