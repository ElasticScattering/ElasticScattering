#include "src/sim/es/cell_grid.h"
#include "settings.h";

__kernel test_get_next_cell_multiple(constant ImpuritySettings *is, constant Orbit *orbits, constant Intersection *starting_intersections, global double2 *next_intersections)
{
	int i = get_global_id(0);

	Intersection *next_intersection;
	bool cell_available = GetNextCell(orbits[i], 0, is, starting_intersections[i], next_intersection);
	next_intersections[i] = next_intersection.position;
}

__kernel test_get_next_cell(constant ImpuritySettings *is, Orbit orbit, Intersection starting_intersection, global double2 *next_intersections)
{
	int i = get_global_id(0);

	Intersection *next_intersection;
	while (1)
	{
		bool cell_available = GetNextCell(orbits[i], 0, is, starting_intersections[i], next_intersection);
		next_intersections[i] = next_intersection.position;
		i++;
	}
}