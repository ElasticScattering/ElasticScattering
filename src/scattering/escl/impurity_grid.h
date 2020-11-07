#pragma once

#ifndef DEVICE_PROGRAM
	#include "src/scattering/escl/v2.h"
#endif

struct Intersection {
	double2 p1, p2;
};

struct CellRange {
	int start, end;
};

double GetAngle(double2 pos, double2 circle, double radius);

Intersection BoundaryIntersects(double2 pos, double2 pos2, double2 circle, double radius, double L); 

void BoxIntersects(double2 low, double L, double2 circle, double circle_radius);

int get_cell_index(const v2 pos, const v2 range, const int cells_per_row);

CellRange GetNextCell();

////////////////////////////

int get_cell_index(const v2 pos, const v2 range, const int cells_per_row);

CellRange get_cell_range(const v2 pos, const v2 range, const int cells_per_row);

////////////////////////


int to_grid(const double x, const double2 range, const int cells_per_row);

v2i to_grid(const double x, const double y, const double2 range, const int cells_per_row);

v2 to_world(const int cell_index, const int cells_per_row, const double2 spawn_range);

bool within_bounds(v2i p, const int cells_per_row);

int to_index(const int x, const int y, const int cells_per_row);

int to_index(const v2i p, const int cells_per_row);

double remap(double x, double s0, double s1, double t0, double t1);
