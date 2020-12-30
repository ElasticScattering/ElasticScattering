#pragma once
#include <doctest.h>
#include "TestMacros.h"
#include "TestUtils.h"

#include "src/scattering/Grid.h"
#include "src/scattering/escl/cell_grid.h"

//Todo: test test_impurities.

int get_cell_index(const double2 pos, const double spawn_region_start, const double spawn_region_size, const int cells_per_row)
{
	return get_index(get_cell(pos, spawn_region_start, spawn_region_size, cells_per_row), cells_per_row);
}

bool ImpurityOverlapsCell(int correct_cell, v2 pos, double impurity_radius, const double spawn_region_start, const double spawn_region_size, int cells_per_row)
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
		auto new_cell = get_cell_index(possible_overlapping_positions[i], spawn_region_start, spawn_region_size, cells_per_row);

		if (correct_cell == new_cell) {
			return true;
		}
	}

	return false;
}

TEST_CASE("Index should be strictly increasing.")
{
	auto grid = Grid(0, 1.0, 0.0, 1e2, 1e-2, 1e2 / 16);
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
	double impurity_radius = 1e-3;
	
	auto grid = Grid(0, 1.0, 0.0, 1e2, impurity_radius, 1e2 / 16.0);
	auto impurity_settings = grid.GetSettings();
	auto impurities = grid.GetImpurities();
	
	double start = impurity_settings.spawn_region_start;
	double size  = impurity_settings.spawn_region_size;
	int cells_per_row = grid.GetCellsPerRow();

	bool in_order = true;
	bool last_one_overlapped = false;
	for (int i = 1; i < impurities.size(); i++)
	{
		auto last_cell_index = get_cell_index(impurities[i - 1], start, size, cells_per_row);
		auto new_cell_index  = get_cell_index(impurities[i], start, size, cells_per_row);

		if (last_cell_index > new_cell_index)
		{
			if (last_one_overlapped)
			{
				last_cell_index = get_cell_index(impurities[i - 2], start, size, cells_per_row);
			}
			if (last_cell_index > new_cell_index) {
				bool overlaps = ImpurityOverlapsCell(last_cell_index, impurities[i], impurity_radius, start, size, cells_per_row);
				if (overlaps) {
					last_one_overlapped = true;
				}
				else {
					in_order = false;
					//printf("Out of order! (%d): !! %d < %d, (%f,%f) < (%f,%f)\n", i, last_cell_index, new_cell_index, impurities[i - 1].x, impurities[i - 1].y, impurities[i].x, impurities[i].y);
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
		auto grid = Grid(0, 1.0, 0.0, 1e2, 1e-3, 4);
		CHECK(grid.GetIndex().size() == 25);
	}

	SUBCASE("One cell should get all impurities.")
	{
		int cells_per_row = 1;
		auto grid = Grid(0, 1.0, 0.0, 1e2, 1e-3, 1e2+1);
		auto index = grid.GetIndex();

		CHECK(index.size() == 1);
		CHECK(index[0] == 100);
	}
}

TEST_CASE("Impurities should be indexed correctly.")
{
	double impurity_radius = 1e-3;
	
	auto grid = Grid(0, 1.0, 0.0, 1e2, impurity_radius, 1e2/16.0);
	auto impurity_settings = grid.GetSettings();

	auto impurities = grid.GetImpurities();
	auto index = grid.GetIndex();

	bool in_order = true;
	for (int i = 0; i < index.size(); i++)
	{
		int impurity_start = (i - 1 < 0) ? 0 : index[i - 1];
		int impurity_end = index[i];

		for (int j = impurity_start; j < impurity_end; j++)
		{
			int actual_cell_index = get_cell_index(impurities[j], impurity_settings.spawn_region_start, impurity_settings.spawn_region_size, impurity_settings.cells_per_row);

			if (i != actual_cell_index && !ImpurityOverlapsCell(i, impurities[j], impurity_radius, impurity_settings.spawn_region_start, impurity_settings.spawn_region_size, impurity_settings.cells_per_row))
			{
				in_order = false;
				//printf("Wrong index! '%d' should be '%d': (%f,%f).\n", i, actual_cell_index, impurities[j].x, impurities[j].y);
				break;
			}
		}
	}

	CHECK(in_order == true);
}



TEST_CASE("Impurities should be indexed diagonally.")
{
	SUBCASE("Impurity in the middle.")
	{
		auto impurities = std::vector<v2>{ { 0.5, 0.5 }, { 0.1, 0.1 }, { 0.9, 0.9 }, { 0.9, 0.8 } };
		auto grid = Grid(impurities, 1, 0, 1e-2, 1);
		REQUIRE_MESSAGE(grid.GetCellsPerRow() == 2, "Niet genoeg cellen om te testen");

		CHECK(grid.GetTotalImpurityCount() == 7);

		/*
		for (int i = 0; i < grid.GetImpurities().size(); i++)
		{
			auto imp = grid.GetImpurities()[i];
			printf("Impurity: %f, %f\n", imp.x, imp.y);
		}
		*/
	}
	
	SUBCASE("Impurity off centered.")
	{
		auto impurities = std::vector<v2>{ { 0.505, 0.4999 }, { 0.1, 0.1 }, { 0.9, 0.9 }, { 0.9, 0.8 } };
		auto grid = Grid(impurities, 1, 0, 1e-2, 1);
		REQUIRE_MESSAGE(grid.GetCellsPerRow() == 2, "Niet genoeg cellen om te testen");

		CHECK(grid.GetTotalImpurityCount() == 7);
	}

	SUBCASE("Impurity off centered low right.")
	{
		auto impurities = std::vector<v2>{ { 0.508, 0.4999 }, { 0.1, 0.1 }, { 0.9, 0.9 }, { 0.9, 0.8 } };
		auto grid = Grid(impurities, 1, 0, 1e-2, 1);
		REQUIRE_MESSAGE(grid.GetCellsPerRow() == 2, "Niet genoeg cellen om te testen");

		CHECK(grid.GetTotalImpurityCount() == 7);
	}

	SUBCASE("Impurity off centered low left.")
	{
		auto impurities = std::vector<v2>{ { 0.492, 0.4999 }, { 0.1, 0.1 }, { 0.9, 0.9 }, { 0.9, 0.8 } };
		auto grid = Grid(impurities, 1, 0, 1e-2, 1);
		REQUIRE_MESSAGE(grid.GetCellsPerRow() == 2, "Niet genoeg cellen om te testen");

		CHECK(grid.GetTotalImpurityCount() == 7);
	}

	SUBCASE("Impurity off centered top right.")
	{
		auto impurities = std::vector<v2>{ { 0.508, 0.501 }, { 0.1, 0.1 }, { 0.9, 0.9 }, { 0.9, 0.8 } };
		auto grid = Grid(impurities, 1, 0, 1e-2, 1);
		REQUIRE_MESSAGE(grid.GetCellsPerRow() == 2, "Niet genoeg cellen om te testen");

		CHECK(grid.GetTotalImpurityCount() == 7);
	}

	SUBCASE("Impurity off centered top left.")
	{
		auto impurities = std::vector<v2>{ { 0.492, 0.501 }, { 0.1, 0.1 }, { 0.9, 0.9 }, { 0.9, 0.8 } };
		auto grid = Grid(impurities, 1, 0, 1e-2, 1);
		REQUIRE_MESSAGE(grid.GetCellsPerRow() == 2, "Niet genoeg cellen om te testen");

		CHECK(grid.GetTotalImpurityCount() == 7);
	}
}