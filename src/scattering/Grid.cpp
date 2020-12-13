#include "Grid.h"

#include "escl/constants.h"

#include <random>
#include <windows.h>
#include <set>

Grid::Grid(int count, int seed, v2 _spawn_range, double impurity_radius, int _cells_per_row)
{
	spawn_range = _spawn_range;
	cells_per_row = _cells_per_row;
	impurity_count = count;

	auto impurities = GenerateImpurities(count, seed);
	GenerateImpurityCells(impurities, impurity_radius);
	ConvertToIndex();
}

Grid::Grid(std::vector<v2> impurities, v2 _spawn_range, double impurity_radius, int _cells_per_row)
{
	spawn_range = _spawn_range;
	cells_per_row = _cells_per_row;
	impurity_count = impurities.size();

	GenerateImpurityCells(impurities, impurity_radius);
	ConvertToIndex();
}


std::vector<v2> Grid::GenerateImpurities(int count, int seed)
{
	std::vector<v2> impurities(count);

	std::uniform_real_distribution<double> unif(spawn_range.x, spawn_range.y);

	std::random_device random_device;
	std::default_random_engine re(seed);

	for (int i = 0; i < count; i++) {
		impurities[i] = { unif(re), unif(re) };
	}

	return impurities;
}

void Grid::GenerateImpurityCells(std::vector<v2> impurities, double impurity_radius)
{
	cells.resize(cells_per_row * cells_per_row);

	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			cells[j * cells_per_row + i] = Cell(i, j);
		}
	}

	for (int i = 0; i < impurities.size(); i++) {
		total_indexed_impurities += add_to_overlapping_cells(cells, impurities[i], impurity_radius);
	}
}

// Move impurities back to single array and build an index.
void Grid::ConvertToIndex()
{
	ordered_impurities.resize(total_indexed_impurities);
	imp_index.resize(cells_per_row * cells_per_row);

	int impurity_counter = 0;
	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			auto imps = cells[j * cells_per_row + i].impurities;

			imp_index[j * cells_per_row + i] = impurity_counter + imps.size();

			for (int k = 0; k < imps.size(); k++) {
				ordered_impurities[impurity_counter] = imps[k];
				impurity_counter++;
			}
		}
	}
}

int Grid::add_to_overlapping_cells(std::vector<Cell>& cells, const v2 pos, const double impurity_radius)
{
	v2 offset = v2(cos(45.0 * PI / 180.0), sin(45.0 * PI / 180.0)) * impurity_radius;

	std::vector<v2i> possible_overlapping_positions = {
		get_cell(pos.x, pos.y					  ),
		get_cell(pos.x + impurity_radius, pos.y	  ),
		get_cell(pos.x - impurity_radius, pos.y	  ),
		get_cell(pos.x, pos.y + impurity_radius	  ),
		get_cell(pos.x, pos.y - impurity_radius	  ),
		get_cell(pos.x + offset.x, pos.y + offset.y),
		get_cell(pos.x - offset.x, pos.y + offset.y),
		get_cell(pos.x + offset.x, pos.y - offset.y),
		get_cell(pos.x - offset.x, pos.y - offset.y)
	};

	std::set<v2i> added_cells;
	for (int i = 0; i < possible_overlapping_positions.size(); i++) {
		auto new_cell = possible_overlapping_positions[i];

		if (added_cells.find(new_cell) == added_cells.end() && within_bounds(new_cell)) {
			added_cells.insert(new_cell);

			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
		}
	}
	
	return added_cells.size();
}

v2i Grid::get_cell(const double x, const double y)
{
	return {
		(int)((x - spawn_range.x) / (spawn_range.y - spawn_range.x) * (double)(cells_per_row)),
		(int)((y - spawn_range.x) / (spawn_range.y - spawn_range.x) * (double)(cells_per_row))
	};
}

bool Grid::within_bounds(const v2i p)
{
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}
