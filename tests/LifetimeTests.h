#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "TestUtils.h"

#include "src/scattering/Grid.h"
#include "src/scattering/escl/lifetime.h"
#include "src/sim_main.h"

void DecideFinalParameters(ScatteringParameters& sp)
{
	if (sp.is_incoherent == 1) sp.tau = HBAR / (KB * sp.temperature);
	sp.default_max_lifetime = 15.0 * sp.tau;

	{
		bool incoherent = (sp.is_incoherent == 1);

		const double incoherent_area = sp.alpha * 2.0;
		sp.integrand_angle_area = incoherent ? incoherent_area : (PI / 2.0 - incoherent_area);
		sp.integrand_step_size = sp.integrand_angle_area / (sp.integrand_steps - 1);

		sp.integrand_start_angle = (incoherent ? -sp.alpha : sp.alpha);
	}
}

TEST_CASE("Lifetime tests")
{
	auto impurities = GetTestImpurities();
	REQUIRE(impurities.size() > 0);

	ScatteringParameters sp;
	sp.alpha = PI / 4.0;

	sp.max_expected_impurities_in_cell = 10;
	sp.impurity_count = impurities.size();
	sp.region_extends = 1e-6;
	sp.region_size = 1e-6;
	sp.impurity_radius = 1e-8;
	sp.impurity_spawn_range = { -sp.region_extends, sp.region_size + sp.region_extends };
	double t = sqrt(sp.impurity_count / (double)sp.max_expected_impurities_in_cell);
	sp.cells_per_row = (int)ceil(sqrt(sp.impurity_count / (double)sp.max_expected_impurities_in_cell));
	sp.cell_size = (sp.impurity_spawn_range.y - sp.impurity_spawn_range.x) / (double)sp.cells_per_row;

	sp.particle_speed = 1.68e5;
	sp.magnetic_field = 10;
	sp.angular_speed = E * sp.magnetic_field / M;

	Grid grid(impurities, sp.impurity_spawn_range, sp.impurity_radius, sp.cells_per_row);

	SUBCASE("Intersect in the second box")
	{
		sp.is_incoherent = 0;
		sp.is_clockwise = 1;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());
		
		CHECK(lt > 0.88);
		CHECK(lt < 0.89);
	}

	SUBCASE("Intersect in the second box, incoherent -> boundtime")
	{
		sp.is_incoherent = 1;
		sp.is_clockwise = 1;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK_ALMOST(lt, (PI / 4));
	}

	SUBCASE("Intersect in first box")
	{
		sp.is_incoherent = 0;
		sp.is_clockwise = 0;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt > 0.078);
		CHECK(lt < 0.079);
	}

	SUBCASE("Boundary limited.")
	{
		sp.is_incoherent = 1;
		sp.is_clockwise = 1;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(7e-7, 2.99e-7), &sp, grid.GetImpurities(), grid.GetIndex());

		double r = abs(lt - PI / 4);
		CHECK(r < 1e-7);
	}

	SUBCASE("Coherent, late intersect.")
	{
		sp.is_incoherent = 0;
		sp.is_clockwise = 1;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(7e-7, 4.3e-7), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt > 3.82); 
		CHECK(lt < 3.83);
	}

	SUBCASE("Repositioned.")
	{
		sp.magnetic_field = 50;
		sp.angular_speed = E * sp.magnetic_field / M;
		sp.is_incoherent = 0;
		sp.is_clockwise = 0;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(6.5e-7, 1.35e-6), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt > 100);
	}

	SUBCASE("Repositioned.")
	{
		sp.magnetic_field = 50;
		sp.angular_speed = E * sp.magnetic_field / M;
		sp.is_incoherent = 0;
		sp.is_clockwise = 0;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;

		double lt = sp.angular_speed * lifetime(0, 0, v2(6.15e-7, 1.35e-6), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt > 2.21);
		CHECK(lt < 2.215);
	}

	SUBCASE("Edge case 1.")
	{
		sp.magnetic_field = 8;
		sp.angular_speed = E * sp.magnetic_field / M;
		sp.is_incoherent = 0;
		sp.is_clockwise = 0;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;
		sp.integrand_step_size = 1.57079633e+00;

		double lt = sp.angular_speed * lifetime(0, 1, v2(1.00000000e-06, 5.0000000e-07), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt > 1);
	}

	SUBCASE("Edge case 2.")
	{
		sp.magnetic_field = 8;
		sp.angular_speed = E * sp.magnetic_field / M;
		sp.is_incoherent = 0;
		sp.is_clockwise = 0;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;
		sp.integrand_step_size = 1.57079633e+00;

		double lt = sp.angular_speed * lifetime(0, 1, v2(1.00000000e-06, 5.0000100e-07), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt < 1);
	}

	SUBCASE("Edge case 3.")
	{
		sp.magnetic_field = 8;
		sp.angular_speed = E * sp.magnetic_field / M;
		sp.is_incoherent = 0;
		sp.is_clockwise = 0;
		DecideFinalParameters(sp);
		sp.integrand_start_angle = 0;
		sp.integrand_step_size = 1.57079633e+00;

		double lt = sp.angular_speed * lifetime(0, 1, v2(1.00000000e-06, 4.9999900e-07), &sp, grid.GetImpurities(), grid.GetIndex());

		CHECK(lt < 1);
	}
}

