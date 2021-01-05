#include "Grid.h"

#include "es/constants.h"
#include "src/utils/ErrorMacros.h"
#include <random>
#include <windows.h>
#include <set>

Grid::Grid(unsigned int seed, double region_size, double region_extends, double density, double _impurity_radius, int target_cell_population)
{
	seed_used = seed;
	impurity_radius = _impurity_radius;
	
	spawn_range = v2(-region_extends, region_size + region_extends);
	double area_length = spawn_range.y - spawn_range.x;
	double imp_count_proposal = area_length * area_length * density;
	CFG_EXIT_CONDITION(imp_count_proposal > 1e7, "Woah!");
	CFG_EXIT_CONDITION(imp_count_proposal < 1, "No impurities");

	unique_impurity_count = max(1, (long)ceil(imp_count_proposal));
	cells_per_row = max((int)round(sqrt(unique_impurity_count / (double)target_cell_population)), 1);
	cell_size = (spawn_range.y - spawn_range.x) / (double)cells_per_row;
	auto impurities = GenerateImpurities(unique_impurity_count, seed);

	GenerateImpurityCells(impurities, impurity_radius);
	ConvertToIndex();
}

Grid::Grid(std::vector<v2> impurities, double region_size, double region_extends, double _impurity_radius, int target_cell_population)
{
	impurity_radius = _impurity_radius;
	
	spawn_range = v2(-region_extends, region_size + region_extends);
	unique_impurity_count = impurities.size();
	cells_per_row = max((int)round(sqrt(unique_impurity_count / (double)target_cell_population)), 1);
	cell_size = (spawn_range.y - spawn_range.x) / (double)cells_per_row;

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
	std::vector<v2i> possible_overlapping_positions = {
		get_cell(pos.x, pos.y),
		get_cell(pos.x + impurity_radius, pos.y),
		get_cell(pos.x - impurity_radius, pos.y),
		get_cell(pos.x, pos.y + impurity_radius),
		get_cell(pos.x, pos.y - impurity_radius),
	};

	std::set<v2i> added_cells;
	for (int i = 0; i < possible_overlapping_positions.size(); i++) {
		auto new_cell = possible_overlapping_positions[i];

		if (added_cells.find(new_cell) == added_cells.end() && within_bounds(new_cell)) {
			added_cells.insert(new_cell);

			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
		}
	}

	int cells_added = added_cells.size();

	v2 low_left  = to_world(possible_overlapping_positions[0]);
	v2 low_right = low_left + v2(cell_size, 0);
	v2 top_left  = low_left + v2(0, cell_size);
	v2 top_right = low_left + v2(cell_size, cell_size);

	double ir2 = impurity_radius * impurity_radius;
	auto cell = get_cell(pos.x, pos.y);

	/*
	std::vector<v2> possible_overlapping_diagonals = {
		low_left,
		low_left + v2(cell_size, 0),
		low_left + v2(0, cell_size),
		low_left + v2(cell_size, cell_size)
	};

	for (int i = 0; i < possible_overlapping_diagonals.size(); i++) {
		auto cell_world_pos = possible_overlapping_diagonals[i];
		double d_squared = pow(cell_world_pos.x - pos.x, 2) + pow(cell_world_pos.y - pos.y, 2);
		if (d_squared < ir2) {
			auto new_cell = get_cell(cell_world_pos.x, cell_world_pos.y);
			if (within_bounds(v2i(new_cell))) {
				cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
				cells_added++;
			}
		}
	}
	*/
	
	double d_topleft_squared = pow(top_left.x - pos.x, 2) + pow(top_left.y - pos.y, 2);
	if (d_topleft_squared < ir2) {
		auto new_cell = v2i(cell.x - 1, cell.y + 1);
		if (within_bounds(v2i(new_cell))) {
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
			cells_added++;
		}
	}

	double d_topright_squared = pow(top_right.x - pos.x, 2) + pow(top_right.y - pos.y, 2);
	if (d_topright_squared < ir2) {
		auto new_cell = v2i(cell.x + 1, cell.y + 1);
		if (within_bounds(v2i(new_cell))) {
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
			cells_added++;
		}
	}

	double d_lowleft_squared = pow(low_left.x - pos.x, 2) + pow(low_left.y - pos.y, 2);
	if (d_lowleft_squared < ir2) {
		auto new_cell = v2i(cell.x - 1, cell.y - 1);
		if (within_bounds(v2i(new_cell))) {
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
			cells_added++;
		}
	}

	double d_lowright_squared = pow(low_right.x - pos.x, 2) + pow(low_right.y - pos.y, 2);
	if (d_lowright_squared < ir2) {
		auto new_cell = v2i(cell.x + 1, cell.y - 1);
		if (within_bounds(v2i(new_cell))) {
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
			cells_added++;
		}
	}

	return cells_added;
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
