#pragma once

#ifndef DEVICE_PROGRAM
	#include "src/scattering/escl/v2.h"
#endif
#include "details.h"

typedef struct Intersection {
	double2 position;
	double incident_angle;

	double dphi;
} Intersection;

typedef struct SideIntersection {
	Intersection i1, i2;
} SideIntersection;

struct CellRange {
	int start, end;
};

double GetAngle(double2 pos, Orbit o);

void FindExitIntersect(Orbit o, double phi, Intersection start_intersect);
bool GetBoundaryIntersects(const double2 p1, const double2 p2, const Orbit o, const double L, SideIntersection* intersection);
bool GetNextCell(const int current_cell, const Orbit orbit, const double2 last_intersection, const int cells_per_row, const double L, const double2 spawn_range, int* next_cell, double2* intersection_point);

int get_cell_index(const v2 pos, const v2 range, const int cells_per_row);

////////////////////////

int to_grid(const double x, const double2 range, const int cells_per_row);

v2i to_grid(const double x, const double y, const double2 range, const int cells_per_row);

v2 to_world(const int cell_index, const int cells_per_row, const double2 spawn_range);

bool within_bounds(v2i p, const int cells_per_row);

int to_index(const int x, const int y, const int cells_per_row);

int to_index(const v2i p, const int cells_per_row);

double remap(double x, double s0, double s1, double t0, double t1);
