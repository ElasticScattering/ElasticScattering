#pragma once

#include "escl/v2.h"
#include "escl/ScatteringParameters.h"

#include <vector>
struct Cell {
	int x;
	int y;

	Cell()
	{
		x = 0;
		y = 0;
	}

	Cell(int _x, int _y)
	{
		x = _x;
		y = _y;
	}

	std::vector<v2> impurities;
};

class ImpurityGridIndex {
private:
	static int add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const v2 spawn_range, const double impurity_radius, const int cells_per_row);
	static v2i to_grid(const double x, double y, const v2 range, const int cells_per_row);
	static bool within_bounds(const v2i p, const int cells_per_row);

public:
	std::vector<v2> impurities;
	std::vector<int> imp_index;

	static ImpurityGridIndex& Generate(int count, int seed, v2 spawn_range, double impurity_radius, int cells_per_row);
};