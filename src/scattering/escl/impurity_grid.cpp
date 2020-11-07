#pragma once

#include "impurity_grid.h"

#ifndef DEVICE_PROGRAM
	#include <cmath>
#endif
#include "constants.h"

double GetAngle(double2 pos, double2 circle, double radius) {
	if (abs(pos.y - circle.y) > EPSILON) {
		double sign = (pos.x < circle.x) ? -1.0 : 1.0;
		double angle = sign * asin((pos.y - circle.y) / radius);
		return fmod(angle, PI2);
	}
	else {
		return (pos.x > circle.x) ? 0 : PI;
	}
}

Intersection BoundaryIntersects(double2 pos, double2 pos2, double2 circle, double radius, double L) {
	double2 u = (pos2 - pos) / L;
	double projection_distance = u.x * (circle.x - pos.x) + u.y * (circle.y - pos.y);
	double2 proj = pos + (u * projection_distance);
	double proj_circle_distance_sq = pow(proj.x - circle.x, 2) + pow(proj.x - circle.x, 2);
	double r2 = (radius * radius);
	Intersection i;
	if (proj_circle_distance_sq >= r2) {
		return i; //should be empty
	}

	double dist = sqrt(r2 - proj_circle_distance_sq);

	double2 intersection_a = proj + u * dist;
	double2 intersection_b = proj - u * dist;

	bool valid_intersect_a = false;
	bool valid_intersect_b = false;

	if (abs(u.x) > abs(u.y)) {
		valid_intersect_a = ((pos.x < intersection_a.x) && (pos2.x > intersection_a.x)) || ((pos.x > intersection_a.x) && (pos2.x < intersection_a.x));
		valid_intersect_b = ((pos.x < intersection_b.x) && (pos2.x > intersection_b.x)) || ((pos.x > intersection_b.x) && (pos2.x < intersection_b.x));
	}
	else {
		valid_intersect_a = ((pos.y < intersection_a.y) && (pos2.y > intersection_a.y)) || ((pos.y > intersection_a.y) && (pos2.y < intersection_a.y));
		valid_intersect_b = ((pos.y < intersection_b.y) && (pos2.y > intersection_b.y)) || ((pos.y > intersection_b.y) && (pos2.y < intersection_b.y));
	}

	if (valid_intersect_a) {
		i.p1 = intersection_a;
	}
	if (valid_intersect_b) {
		i.p2 = intersection_b;
	}
}

void BoxIntersects(double2 low, double L, double2 circle, double circle_radius) {
	Intersection right = BoundaryIntersects(low, low + double2(L, 0), circle, circle_radius, L);
	Intersection up = BoundaryIntersects(low, low + double2(0, L), circle, circle_radius, L);
	Intersection down = BoundaryIntersects(low + double2(L, 0), low + double2(L, L), circle, circle_radius, L);
	Intersection left = BoundaryIntersects(low + double2(0, L), low + double2(L, L), circle, circle_radius, L);

	//Intersects = <filter out the Empty>
	//Intersects3 = AddAngleToIntersects(Intersects, xc, yc, rc)
	//Return Intersects3
}

CellRange GetNextCell(int current_cell, int cells_per_row, double2 spawn_range, double2 orbit, double orbit_radius)
{
	v2 low_left = to_world(current_cell, cells_per_row, spawn_range);

	double L = (spawn_range.y - spawn_range.x) / (double)cells_per_row;
	auto intersects = BoxIntersects(low_left, L, );
}

////////////////////////////

int get_cell_index(const v2 pos, const v2 range, const int cells_per_row)
{
	return to_index(to_grid(pos.x, pos.y, range, cells_per_row), cells_per_row);
}

CellRange get_cell_range(const v2 pos, const v2 range, const int cells_per_row)
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