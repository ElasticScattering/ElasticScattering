#pragma once

#include "datastructures/v2.h"
#include "datastructures/ScatteringParameters.h"

#include <vector>

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

struct ImpuritySettings {
	double region_size;
	double region_extends;
	double cells_per_row;

	double total_impurity_count;

	double impurity_radius;

	unsigned int seed;
};

class ImpurityGrid {
	double impurity_radius;
	v2 range;
	
	int cells_per_row;
	
	std::vector<Cell> cells;

	int to_grid(const double x) const;
	v2i to_grid(const double x, double y) const;
	v2i to_grid(const v2 p) const;
	int to_index(const v2i p) const;
	int to_index(const int x, const int y) const;

	//bool within_bounds(int x, int og_x);
	//bool within_bounds(v2i p, v2i og_p);

	void add_to_overlapping_cells(const v2 pos);

public:
	std::vector<v2> impurities;
	std::vector<int> imp_index;

	ImpurityGrid(ImpuritySettings settings);
};