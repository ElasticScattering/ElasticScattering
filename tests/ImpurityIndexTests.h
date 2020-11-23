#pragma once
#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/ImpurityIndex.h"
#include "src/scattering/escl/impurity_grid.h"

TEST_CASE("Index be strictly increasing")
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
			printf("Out of order! %d should not be less than %d, maar is %d", index[i], index[i-1], index[i] < index[i - 1]);
			break;
		}
	}

	CHECK(in_order);
}

TEST_CASE("Impurities should be indexed in order.")
{
	v2 range = { 0.0, 1.0 };
	int cells_per_row = 4;
	auto grid = ImpurityIndex(100, 0, range, 1e-2, cells_per_row);
	auto impurities = grid.GetImpurities();
	
	bool in_order = true;
	for (int i = 1; i < impurities.size(); i++)
	{
		auto last_cell_index = get_cell_index(impurities[i - 1], range, cells_per_row);
		auto new_cell_index = get_cell_index(impurities[i], range, cells_per_row);

		if (!(last_cell_index <= new_cell_index))
		{
			in_order = false;
			printf("Out of order! %d should not be greater than %d", last_cell_index, new_cell_index);
			break;
		}
	}

	CHECK(in_order);
}