#include "ImpurityGridIndex.h"

#include "escl/constants.h"

#include <random>
#include <windows.h>
#include <set>

ImpurityGridIndex& ImpurityGridIndex::Generate(int count, int seed, v2 spawn_range, double impurity_radius, int cells_per_row)
{
	ImpurityGridIndex grid;

	std::vector<Cell> cells(cells_per_row * cells_per_row);

	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			cells[j * cells_per_row + i] = Cell(i, j);
		}
	}

	int total_impurity_entities = 0;
	// Generate impurities and put them into a cell.
	{
		std::uniform_real_distribution<double> unif(spawn_range.x, spawn_range.y);

		std::random_device random_device;
		std::default_random_engine re(seed);

		for (int i = 0; i < count; i++) {
			v2 pos = { unif(re), unif(re) };

			total_impurity_entities += add_to_overlapping_cells(cells, pos, spawn_range, impurity_radius, cells_per_row);
		}
	}

	// Move impurities to single array and build an index.
	grid.imp_index.resize(max(cells_per_row * cells_per_row, 2));
	grid.impurities.resize(total_impurity_entities);

	int impurity_counter = 0;
	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			auto imps = cells[j * cells_per_row + i].impurities;
			for (int k = 0; k < imps.size(); k++) {
				grid.impurities[impurity_counter] = imps[k];
				impurity_counter++;
			}

			grid.imp_index[j * cells_per_row + i] = imps.size() + impurity_counter;
		}
	}

	return grid;
}

int ImpurityGridIndex::add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const v2 spawn_range, const double impurity_radius, const int cells_per_row)
{
	v2 offset = v2(cos(45.0 * PI / 180.0), sin(45.0 * PI / 180.0)) * impurity_radius;

	std::vector<v2i> possible_overlapping_positions = {
		to_grid(pos.x + impurity_radius, pos.y	  , spawn_range, cells_per_row),
		to_grid(pos.x - impurity_radius, pos.y	  , spawn_range, cells_per_row),
		to_grid(pos.x, pos.y + impurity_radius	  , spawn_range, cells_per_row),
		to_grid(pos.x, pos.y - impurity_radius	  , spawn_range, cells_per_row),
		to_grid(pos.x + offset.x, pos.y + offset.y, spawn_range, cells_per_row),
		to_grid(pos.x - offset.x, pos.y + offset.y, spawn_range, cells_per_row),
		to_grid(pos.x + offset.x, pos.y - offset.y, spawn_range, cells_per_row),
		to_grid(pos.x - offset.x, pos.y - offset.y, spawn_range, cells_per_row)
	};

	// Add main cell.
	v2i cell = to_grid(pos.x, pos.y, spawn_range, cells_per_row);
	cells[cell.y * cells_per_row + cell.x].impurities.push_back(pos);

	std::set<v2i> added_cells = { cell };
	for (int i = 0; i < possible_overlapping_positions.size(); i++) {
		auto new_cell = possible_overlapping_positions[i];

		if (cell != new_cell && added_cells.find(new_cell) == added_cells.end() && within_bounds(new_cell, cells_per_row)) {
			added_cells.insert(new_cell);
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
		}
	}
	
	//printf("Main cell: (%d,%d), total: %d\n", cell.x, cell.y, added_cells.size());

	return added_cells.size();
}

v2i ImpurityGridIndex::to_grid(const double x, const double y, const v2 range, const int cells_per_row)
{
	return {
		(int)((x - range.x) / (range.y - range.x) * (cells_per_row)),
		(int)((y - range.x) / (range.y - range.x) * (cells_per_row))
	};
}

bool ImpurityGridIndex::within_bounds(const v2i p, const int cells_per_row)
{
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}
