#pragma once

#ifndef DEVICE_PROGRAM
#include "src/datastructures/v2.h"
#endif

void BoundaryIntersects(double2 pos, double2 pos2, double2 circle, double radius, double L) {
	double2 u = (pos2 - pos) / L;
	double projection_distance = u.x * (circle.x - pos.x) + u.y * (circle.y - pos.y);
	double2 proj = pos + (u * projection_distance);
	double proj_circle_distance_sq = pow(proj.x - circle.x, 2) + pow(proj.x - circle.x, 2);
	double r2 = (radius * radius);
	if (proj_circle_distance_sq >= r2) {
		return; //empty
	}

	double dist = sqrt(r2 - proj_circle_distance_sq);

	double2 intersection_a = proj + u * dist;
	double2 intersection_b = proj + u * dist;

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

	}
	if (valid_intersect_b) {

	}
}
	

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

bool within_bounds(v2i p, const int cells_per_row) {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}

int to_index(const int x, const int y, const int cells_per_row) {
	return y * cells_per_row + x;
}

int to_index(const v2i p, const int cells_per_row) {
	to_index(p.x, p.y, cells_per_row);
}

double remap(double x, double s0, double s1, double t0, double t1)
{
	return t0 + (x - s0) / (s1 - s0) * (t1 - t0);
}