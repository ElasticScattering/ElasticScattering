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

class ImpurityIndex {
private:
	std::vector<Cell> cells;
	std::vector<v2> impurities;
	std::vector<int> imp_index;

	int total_indexed_impurities;
	double cells_per_row;
	v2 spawn_range;

	int add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const double impurity_radius);
	v2i get_cell(const double x, double y);
	bool within_bounds(const v2i p);

	void GenerateImpurityCells(int count, int seed, double impurity_radius);
	void ConvertToIndex();

public:
	const std::vector<v2>& GetImpurities() const { return impurities; };
	const std::vector<int>& GetIndex() const { return imp_index; };

	ImpurityIndex(int count, int seed, v2 spawn_range, double impurity_radius, int cells_per_row);
};