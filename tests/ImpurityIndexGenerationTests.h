#pragma once
#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/ImpurityIndex.h"
#include "src/scattering/escl/cell_grid.h"

bool ImpurityOverlapsCell(int correct_cell, v2 pos, double impurity_radius, double2 range, int cells_per_row)
{
	v2 offset = v2(cos(45.0 * PI / 180.0), sin(45.0 * PI / 180.0)) * impurity_radius;

	std::vector<v2> possible_overlapping_positions = {
		v2(pos.x + impurity_radius, pos.y),
		v2(pos.x - impurity_radius, pos.y),
		v2(pos.x, pos.y + impurity_radius),
		v2(pos.x, pos.y - impurity_radius),
		v2(pos.x + offset.x, pos.y + offset.y),
		v2(pos.x - offset.x, pos.y + offset.y),
		v2(pos.x + offset.x, pos.y - offset.y),
		v2(pos.x - offset.x, pos.y - offset.y)
	};

	for (int i = 0; i < possible_overlapping_positions.size(); i++) {
		auto new_cell = get_cell_index(possible_overlapping_positions[i], range, cells_per_row);

		if (correct_cell == new_cell) {
			//printf("(%f,%f) was %d and should be %d.\n", possible_overlapping_positions[i].x, possible_overlapping_positions[i].y, new_cell, correct_cell);
			return true;
		}
		//printf("(%f,%f) was %d but should be %d.\n", possible_overlapping_positions[i].x, possible_overlapping_positions[i].y, new_cell, correct_cell);
	}

	return false;
}

TEST_CASE("Index should be strictly increasing.")
{
	v2 range = { 0.0, 1.0 };
	int cells_per_row = 4;
	auto grid = ImpurityIndex(100, 0, range, 1e-2, cells_per_row);
	auto index = grid.GetIndex();

	bool in_order = true;
	for (int i = 1; i < index.size(); i++)
	{
		if (index[i] < index[i-1])
		{
			in_order = false;
			break;
		}
	}

	CHECK(in_order);
}

TEST_CASE("Impurities should be indexed in order.")
{
	v2 range = { 0.0, 1.0 };
	int cells_per_row = 4;
	double impurity_radius = 1e-3;
	auto grid = ImpurityIndex(100, 0, range, impurity_radius, cells_per_row);
	auto impurities = grid.GetImpurities();
	
	bool in_order = true;
	bool last_one_overlapped = false;
	for (int i = 1; i < impurities.size(); i++)
	{
		auto last_cell_index = get_cell_index(impurities[i - 1], range, cells_per_row);
		auto new_cell_index  = get_cell_index(impurities[i], range, cells_per_row);

		if (last_cell_index > new_cell_index)
		{
			if (last_one_overlapped)
			{
				last_cell_index = get_cell_index(impurities[i - 2], range, cells_per_row);
			}
			if (last_cell_index > new_cell_index) {
				bool overlaps = ImpurityOverlapsCell(last_cell_index, impurities[i], impurity_radius, range, cells_per_row);
				if (overlaps) {
					last_one_overlapped = true;
				}
				else {
					in_order = false;
					printf("Out of order! (%d): !! %d < %d, (%f,%f) < (%f,%f)", i, last_cell_index, new_cell_index, impurities[i - 1].x, impurities[i - 1].y, impurities[i].x, impurities[i].y);
					break;
				}
			}
		}
	}

	CHECK(in_order == true);
}

TEST_CASE("Index should contain all cells.")
{
	SUBCASE("Regular number of cells per row should get squared amount.")
	{
		int cells_per_row = 5;
		auto grid = ImpurityIndex(100, 0, { 0.0, 1.0 }, 1e-3, cells_per_row);

		CHECK(grid.GetIndex().size() == 25);
	}

	SUBCASE("One cell should get all impurities.")
	{
		int cells_per_row = 1;
		auto grid = ImpurityIndex(100, 0, { 0.0, 1.0 }, 1e-3, cells_per_row);
		auto index = grid.GetIndex();

		CHECK(index.size() == 1);
		CHECK(index[0] == 100);
	}
}

TEST_CASE("Impurities should be indexed correctly.")
{
	v2 range = { 0.0, 1.0 };
	int cells_per_row = 5;
	double impurity_radius = 1e-3;
	auto grid = ImpurityIndex(100, 0, range, impurity_radius, cells_per_row);
	auto impurities = grid.GetImpurities();
	auto index = grid.GetIndex();

	bool in_order = true;
	for (int i = 0; i < index.size(); i++)
	{
		int impurity_start = (i - 1 < 0) ? 0 : index[i - 1];
		int impurity_end = index[i];

		for (int j = impurity_start; j < impurity_end; j++)
		{
			int actual_cell_index = get_cell_index(impurities[j], range, cells_per_row);

			if (i != actual_cell_index && !ImpurityOverlapsCell(i, impurities[j], impurity_radius, range, cells_per_row))
			{
				in_order = false;
				printf("Wrong index! '%d' should be '%d': (%f,%f).", i, actual_cell_index, impurities[j].x, impurities[j].y);
				break;
			}
		}
	}

	CHECK(in_order == true);
}
