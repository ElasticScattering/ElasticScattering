#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/escl/lifetime.h"

TEST_CASE("Lifetime tests")
{
	/*
	region, extends = 1e-6, 1e-6
    particle = np.array([7e-7, 2.99e-7, 0])
    vf = 1.68e5
    wc = E * 10 / (5 * 9.109e-31)
    imps = get_test_impurities()
    imp3d, nr1D = grid_divide(imps, 1e-8, region, extends, 10)
    alpha = np.pi / 4
    imp_rad = 1e-8

    # Intersect in the second box
    t = lifetime_full(particle, vf, wc,
                      region, extends, nr1D, imp3d, imp_rad,
                      alpha, True, True)
    radians = wc * t
    assert(0.88 < radians < 0.89)

	*/
	ScatteringParameters sp;
	sp.impurity_radius = 1e-8;
	sp.alpha = PI / 4.0;
	sp.region_extends = 1e-6;
	sp.region_size = 1e-6;

	double wc = E * 2.5 / M0;

	CellInde

	SUBCASE("Intersect in the second box")
	{
		double lt = wc * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp,  );

		CHECK(lt > 0.88 && lt < 0.89);
	}
}

