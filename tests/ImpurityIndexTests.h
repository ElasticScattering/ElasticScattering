#pragma once
#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/ImpurityIndex.h"


TEST_CASE("Index should be ordered correctly")
{
	auto grid = ImpurityIndex(100, 0, { 0.0, 1.0 }, 1e-2, 4);
	auto impurities = grid.GetImpurities();

	bool in_order = true;
	for (int i = 1; i < impurities.size(); i++)
	{
		if (!(impurities[i - 1] < impurities[i]))
		{
			in_order = false;
			break;
		}
	}

	CHECK(in_order);
}