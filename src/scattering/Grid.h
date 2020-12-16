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

class Grid {
private:
	std::vector<Cell> cells;
	std::vector<v2> ordered_impurities;
	std::vector<int> imp_index;

	int unique_impurity_count;
	int total_indexed_impurities;
	int cells_per_row;
	v2 spawn_range;

	int add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const double impurity_radius);
	v2i get_cell(const double x, double y);
	bool within_bounds(const v2i p);

	std::vector<v2> GenerateImpurities(int count, int seed);
	void GenerateImpurityCells(std::vector<v2> impurities, double impurity_radius);
	void ConvertToIndex();

public:
	const std::vector<v2>& GetImpurities() const { return ordered_impurities; };
	const std::vector<int>& GetIndex() const { return imp_index; };

	Grid(int count, int seed, v2 spawn_range, double impurity_radius, int cells_per_row);
	Grid(std::vector<v2> impurities, v2 _spawn_range, double impurity_radius, int _cells_per_row);
};