#include "ImpurityIndex.h"

#include "escl/constants.h"

#include <random>
#include <windows.h>
#include <set>

ImpurityIndex::ImpurityIndex(int count, int seed, v2 _spawn_range, double impurity_radius, int _cells_per_row)
{
	spawn_range = _spawn_range;
	cells_per_row = _cells_per_row;

	GenerateImpurityCells(count, seed, impurity_radius);
	ConvertToIndex();
}

void ImpurityIndex::GenerateImpurityCells(int count, int seed, double impurity_radius)
{
	cells.resize(cells_per_row * cells_per_row);

	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			cells[j * cells_per_row + i] = Cell(i, j);
		}
	}

	// Generate impurities and put them into a cell.
	{
		std::uniform_real_distribution<double> unif(spawn_range.x, spawn_range.y);

		std::random_device random_device;
		std::default_random_engine re(seed);

		for (int i = 0; i < count; i++) {
			v2 pos = { unif(re), unif(re) };
			total_indexed_impurities += add_to_overlapping_cells(cells, pos, impurity_radius);
		}
	}
}

// Move impurities to single array and build an index.
void ImpurityIndex::ConvertToIndex()
{
	impurities.resize(total_indexed_impurities);
	imp_index.resize(max(cells_per_row * cells_per_row, 2));

	int impurity_counter = 0;
	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			auto imps = cells[j * cells_per_row + i].impurities;
			for (int k = 0; k < imps.size(); k++) {
				impurities[impurity_counter] = imps[k];
				impurity_counter++;
			}

			imp_index[j * cells_per_row + i] = imps.size() + impurity_counter;
		}
	}
}

int ImpurityIndex::add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const double impurity_radius)
{
	//@Todo: check of dit idd diagonale cellen is.
	v2 offset = v2(cos(45.0 * PI / 180.0), sin(45.0 * PI / 180.0)) * impurity_radius;

	std::vector<v2i> possible_overlapping_positions = {
		get_cell(pos.x + impurity_radius, pos.y	  ),
		get_cell(pos.x - impurity_radius, pos.y	  ),
		get_cell(pos.x, pos.y + impurity_radius	  ),
		get_cell(pos.x, pos.y - impurity_radius	  ),
		get_cell(pos.x + offset.x, pos.y + offset.y),
		get_cell(pos.x - offset.x, pos.y + offset.y),
		get_cell(pos.x + offset.x, pos.y - offset.y),
		get_cell(pos.x - offset.x, pos.y - offset.y)
	};

	// Add main cell.
	v2i cell = get_cell(pos.x, pos.y);
	cells[cell.y * cells_per_row + cell.x].impurities.push_back(pos);

	std::set<v2i> added_cells = { cell };
	for (int i = 0; i < possible_overlapping_positions.size(); i++) {
		auto new_cell = possible_overlapping_positions[i];

		if (cell != new_cell && added_cells.find(new_cell) == added_cells.end() && within_bounds(new_cell)) {
			added_cells.insert(new_cell);
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
		}
	}
	
	//printf("Main cell: (%d,%d), total: %d\n", cell.x, cell.y, added_cells.size());
	return added_cells.size();
}

v2i ImpurityIndex::get_cell(const double x, const double y)
{
	return {
		(int)((x - spawn_range.x) / (spawn_range.y - spawn_range.x) * (cells_per_row)),
		(int)((y - spawn_range.x) / (spawn_range.y - spawn_range.x) * (cells_per_row))
	};
}

bool ImpurityIndex::within_bounds(const v2i p)
{
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}
