#include "ImpurityGrid.h"

#include "escl/constants.h"

#include <random>
#include <windows.h>

void ImpurityGrid::Generate(ScatteringParameters& sp) {
	cells_per_row = sp.cells_per_row;

	cells.resize(cells_per_row * cells_per_row);

	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			Cell c(i, j);
			cells[j * cells_per_row + i] = c;
		}
	}

	// Generate impurities and put them into a cell.
	{
		std::uniform_real_distribution<double> unif(sp.impurity_spawn_range.x, sp.impurity_spawn_range.y);

		std::random_device random_device;
		std::default_random_engine re(sp.impurity_seed);

		for (int i = 0; i < sp.impurity_count; i++) {
			v2 pos = { unif(re), unif(re) };

			v2i cell = to_grid(pos);
			cells[cell.y * cells_per_row + cell.x].impurities.push_back(pos);

			add_to_overlapping_cells(pos);
		}
	}

	// Move impurities to single array and build an index.
	imp_index.resize(max(cells_per_row * cells_per_row, 2));

	int impurity_counter = 0;
	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			auto imps = cells[j * cells_per_row + i].impurities;
 			for (int k = 0; k < imps.size(); k++) {
				impurities[impurity_counter++] = imps[k];
			}

			imp_index[j * cells_per_row + i] = imps.size() + impurity_counter;
		}
	}
}

void ImpurityGrid::add_to_overlapping_cells(const v2 pos)
{
	v2i cell = to_grid(pos);
	cells[cell.y * cells_per_row + cell.x].impurities.push_back(pos);

	v2 offset = v2(cos(45.0 * PI / 180.0), sin(45.0 * PI / 180.0)) * impurity_radius;

	std::vector<v2i> possible_overlapping_positions = {
		to_grid(pos.x + impurity_radius, pos.y),
		to_grid(pos.x - impurity_radius, pos.y),
		to_grid(pos.x, pos.y + impurity_radius),
		to_grid(pos.x, pos.y - impurity_radius),
		to_grid(pos.x + offset.x, pos.y + offset.y),
		to_grid(pos.x - offset.x, pos.y + offset.y),
		to_grid(pos.x + offset.x, pos.y - offset.y),
		to_grid(pos.x - offset.x, pos.y - offset.y)
	};

	for (int i = 0; i < possible_overlapping_positions.size(); i++) {
		auto new_cell = possible_overlapping_positions[i];

		if (cell != new_cell && within_bounds(new_cell)) {
			cells[new_cell.y * cells_per_row + new_cell.x].impurities.push_back(pos);
		}
	}
}

v2i ImpurityGrid::to_grid(const v2 p) const
{
	return to_grid(p.x, p.y);
}

v2i ImpurityGrid::to_grid(const double x, const double y) const
{
	return {
		(int)((x - range.x) / (range.y - range.x) * (cells_per_row)),
		(int)((y - range.x) / (range.y - range.x) * (cells_per_row))
	};
}

bool ImpurityGrid::within_bounds(v2i p) const {
	return (p.x >= 0 && p.x < cells_per_row) && (p.y >= 0 && p.y < cells_per_row);
}
