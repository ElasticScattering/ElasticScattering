#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "TestUtils.h"

#include "src/scattering/Grid.h"
#include "src/scattering/escl/lifetime.h"
#include "src/sim_main.h"

TEST_CASE("Lifetime tests")
{
	ScatteringParameters sp;
	sp.alpha = PI / 4.0;
	sp.region_extends = 1e-6;
	sp.region_size = 1e-6;
    sp.impurity_radius = 1e-8;
    sp.max_expected_impurities_in_cell = 10;
    sp.particle_speed = 1.67834e5;

    CompleteSimulationParameters(sp);

    double wc = E * 10 / M;

    auto impurities = GetTestImpurities();
    REQUIRE(impurities.size() > 0);

    Grid grid(impurities, sp.impurity_spawn_range, sp.impurity_radius, sp.cells_per_row);

	SUBCASE("Intersect in the second box")
	{
        sp.is_incoherent = 1;
        sp.is_clockwise = 1;
        double lt = wc * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());

        printf("LT: %f, %e\n", lt, lt);
        CHECK(lt > 0.88);
        CHECK(lt < 0.89);
	}

    SUBCASE("Intersect in first box")
    {
        sp.is_incoherent = 1;
        sp.is_clockwise = 0;
        double lt = wc * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());

        printf("LT: %f, %e\n", lt, lt);
        CHECK(lt > 0.078);
        CHECK(lt < 0.079);
    }

    SUBCASE("3rd test")
    {
        sp.is_incoherent = 0;
        sp.is_clockwise = 1;
        double lt = wc * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());

        printf("LT: %f, %e\n", lt, lt);
        double r = abs(lt - PI / 4);
        CHECK(r < 1e-7);
    }
}

