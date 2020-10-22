#pragma once

#include "datastructures/v2.h"
#include "datastructures/ScatteringParameters.h"

#include <vector>

class ImpurityGrid {
	struct Cell {
		int x;
		int y;

		Cell(int _x, int _y)
		{
			x = _x;
			y = _y;
		}

		std::vector<v2> impurities;
	};

	int cells_per_row;
	
	double impurity_radius;
	v2 range;
	
	std::vector<Cell> cells;

	int to_grid(const double x) const;
	v2i to_grid(const double x, double y) const;
	v2i to_grid(const v2 p) const;
	int to_index(const v2i p) const;
	int to_index(const int x, const int y) const;

	bool within_bounds(v2i p) const;

	void add_to_overlapping_cells(const v2 pos);

public:
	std::vector<v2> impurities;
	std::vector<int> imp_index;

	void Generate(ScatteringParameters& settings);
};