/* -------------------------------------------------------------------------
	This code is part of ElasticScattering.

	Copyright(C) 2022 Elastic Scattering developers

	This program is free software : you can redistribute it and /or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#pragma once

#include "es/v2.h"
#include "es/settings.h"

#include <vector>
struct Cell {
	int x;
	int y;

	Cell(int _x, int _y)
	{
		x = _x;
		y = _y;
	}

	Cell() : Cell(0, 0)
	{}

	std::vector<v2> impurities;
};

// Creates impurities distributed over a field, then builds an index that makes it possible 
// to lookup nearby impurities at a position in the field.
class Grid {
private:
	std::vector<Cell> cells;
	std::vector<v2> ordered_impurities;
	std::vector<int> imp_index;

	long unique_impurity_count = 0;
	int total_indexed_impurities = 0;
	
	v2 spawn_range;
	int cells_per_row = 0;
	double cell_size = 0;

	double impurity_radius = 0;
	unsigned int seed_used = 0;

	int add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const double impurity_radius);
	v2i get_cell(const double x, double y);
	bool within_bounds(const v2i p);

	v2 to_world(const v2i current_cell)
	{
		return v2(spawn_range.x + current_cell.x * cell_size, spawn_range.x + current_cell.y * cell_size);
	}

	std::vector<v2> GenerateImpurities(int count, int seed);
	void GenerateImpurityCells(std::vector<v2> impurities, double impurity_radius);
	void ConvertToIndex();

public:
	const std::vector<v2>& GetImpurities() const { return ordered_impurities; };
	const std::vector<int>& GetIndex() const { return imp_index; };

	int GetUniqueImpurityCount() const { return unique_impurity_count; };
	int GetTotalImpurityCount() const { return total_indexed_impurities; };
	int GetCellsPerRow() const { return cells_per_row; };
	int GetSeed() const { return seed_used; };

	ImpuritySettings GetSettings() const
	{
		ImpuritySettings is;
		is.impurity_radius    = impurity_radius;
		is.spawn_region_start = spawn_range.x;
		is.spawn_region_size  = spawn_range.y - spawn_range.x;
		is.cells_per_row      = cells_per_row;
		is.cell_size          = is.spawn_region_size / (double)cells_per_row;
		return is;
	}

	Grid(unsigned int seed, double region_size, double region_extends, double density, double _impurity_radius, int target_impurity_count_per_cell);
	Grid(std::vector<v2> impurities, double region_size, double region_extends, double impurity_radius, int target_impurity_count_per_cell);
};