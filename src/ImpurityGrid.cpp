#include "ImpurityGrid.h"
#include <random>
#include "src/escl/constants.h"

ImpurityGrid::ImpurityGrid(ImpuritySettings settings) {
	cells_per_row= settings.cells_per_row;
	cells.resize(cells_per_row* cells_per_row);

	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			Cell c(i, j);
			c.impurities.clear();
		}
	}

	v2 range = { -settings.region_extends, settings.region_size + settings.region_extends };
	double total_spawn_size = range.y - range.x;

	{
		std::uniform_real_distribution<double> unif(range.x, range.y);

		std::random_device random_device;
		std::default_random_engine re(settings.seed);

		for (int i = 0; i < settings.total_impurity_count; i++) {
			v2 pos = { unif(re), unif(re) };

			v2i cell = to_grid(pos);
			cells[to_index(cell)].impurities.push_back(pos);

			add_to_overlapping_cells(pos);
		}
	}

	imp_index.resize(cells_per_row* cells_per_row);

	int z = 0;
	for (int j = 0; j < cells_per_row; j++) {
		for (int i = 0; i < cells_per_row; i++) {
			auto imps = cells[j * cells_per_row+ i].impurities;
			imp_index[j * cells_per_row+ i] = imps.size();

			for (int k = 0; k < imps.size(); k++)
			{
				impurities[z++] = imps[k];
			}
		}
	}
}

int ImpurityGrid::to_grid(const double x) const
{
	return (int)((x - range.x) / (range.y - range.x) * (cells_per_row));
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

bool within_bounds(int x, int og_x, int bound) {
	return x != og_x && x >= 0 && x < bound;
}

bool within_bounds(v2i p, v2i og_p, int bound) {
	return !(p.x == og_p.x && p.y == og_p.y) && p.x >= 0 && p.x < bound&& p.y >= 0 && p.y < bound;
}

int ImpurityGrid::to_index(const v2i p) const {
	to_index(p.x, p.y);
}

int ImpurityGrid::to_index(const int x, const int y) const {
	return y * cells_per_row + x;
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
		auto p = possible_overlapping_positions[i];

		if (within_bounds(p, cell, cells_per_row)) {
			cells[to_index(p)].impurities.push_back(pos);
		}
	}
}
