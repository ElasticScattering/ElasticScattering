#pragma once

#ifndef DEVICE_PROGRAM
	#include "src/scattering/escl/v2.h"
#endif
#include "details.h"

typedef struct Intersection {
	double2 position;
	double incident_angle;
	double dphi;
	int entered_cell;
} Intersection;

struct CellRange {
	int start, end;
};


double GetAngle(double2 pos, Orbit o);

//bool GetBoundaryIntersects(const double2 p1, const double2 p2, const Orbit o, const double L, SideIntersection* intersection);
bool GetFirstBoundaryIntersect(const double2 p1, const double2 p2, const Orbit o, const double L, const double start_phi, Intersection* intersection);

bool GetNextCell(const Orbit orbit, const int current_cell, const double2 current_cell_lowleft, const Intersection last_intersection, const double L, const int cells_per_row, int* next_cell, Intersection* next_intersection);

int get_cell_index(const v2 pos, const v2 range, const int cells_per_row);

////////////////////////

int to_grid(const double x, const double2 range, const int cells_per_row);

v2i to_grid(const double x, const double y, const double2 range, const int cells_per_row);

v2 to_world(const int cell_index, const int cells_per_row, const double2 spawn_range);

bool within_bounds(v2i p, const int cells_per_row);

int to_index(const int x, const int y, const int cells_per_row);

int to_index(const v2i p, const int cells_per_row);

double remap(double x, double s0, double s1, double t0, double t1);
