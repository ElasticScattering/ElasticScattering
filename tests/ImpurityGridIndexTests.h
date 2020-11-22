#pragma once
#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/ImpurityGridIndex.h"


TEST_CASE("Index should be ordered correctly")
{
	auto grid = ImpurityGridIndex::Generate(100, 0, { 0.0, 1.0 }, 1e-2, 4);

	bool in_order = true;
	for (int i = 1; i < grid.impurities.size(); i++)
	{
		if (!(grid.impurities[i - 1] < grid.impurities[i]))
		{
			in_order = false;
			break;
		}
	}

	CHECK(in_order);
}