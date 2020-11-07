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

class ImpurityGrid {
	int cells_per_row;
	double impurity_radius;
	v2 range;
	
	std::vector<Cell> cells;

	v2i to_grid(const double x, double y) const;
	v2i to_grid(const v2 p) const;

	bool within_bounds(v2i p) const;

	void add_to_overlapping_cells(const v2 pos);

public:
	std::vector<v2> impurities;
	std::vector<int> imp_index;

	void Generate(ScatteringParameters& settings);
};