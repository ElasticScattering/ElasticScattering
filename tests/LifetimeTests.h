#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/escl/lifetime.h"

TEST_CASE("Lifetime tests")
{

	ScatteringParameters sp;
	sp.impurity_radius = 1e-8;
	sp.alpha = PI / 4.0;
	sp.region_extends = 1e-6;
	sp.region_size = 1e-6;

	double wc = E * 2.5 / M0;

	SUBCASE("Intersect in the second box")
	{
		double lt = wc * lifetime(0, 0, ..., &sp, );

		CHECK(lt > 0.88 && lt < 0.89);
	}
}

