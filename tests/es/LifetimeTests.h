#pragma once

#include <doctest.h>
#include "tests/TestMacros.h"
#include "tests/TestUtils.h"

#include "src/sim/Grid.h"
#include "src/sim/es/lifetime.h"

TEST_CASE("Lifetime tests")
{
	auto impurities = GetTestImpurities("tests/data/test_impurities_lt.dat");
	REQUIRE(impurities.size() > 0);

	auto grid = Grid(impurities, 1e-6, 1e-6, 1e-8, 10);
	auto impurity_settings = grid.GetSettings();

	ParticleMetrics pm;

	ParticleSettings ps;
	ps.alpha          = PI / 4.0;
	ps.particle_speed = 1.68e5;
	ps.phi_start      = 0;
	ps.phi_step_size  = 0;
	ps.angular_speed  = E * 10 / M;

	auto ps50 = ps;
	ps50.angular_speed = E * 50 / M;

	auto ps8 = ps;
	ps8.angular_speed = E * 8 / M;
	ps8.phi_step_size = 1.57079633e+00;

	SUBCASE("Intersect in the second box")
	{
		ps.is_clockwise = 1;
		ps.is_coherent  = 1;

		auto particle = CreateParticle(0, 0, v2(7e-7, 2.99e-7), &ps);
		double lt = ps.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);
		
		CHECK(lt > 0.88);
		CHECK(lt < 0.89);
	}

	SUBCASE("Intersect in the second box, incoherent -> boundtime")
	{
		ps.is_clockwise = 1;
		ps.is_coherent  = 0;

		auto particle = CreateParticle(0, 0, v2(7e-7, 2.99e-7), &ps);
		double lt = ps.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK_ALMOST(lt, (PI / 4));
	}

	SUBCASE("Intersect in first box")
	{
		ps.is_clockwise = 0;
		ps.is_coherent  = 1;

		auto particle = CreateParticle(0, 0, v2(7e-7, 2.99e-7), &ps);
		double lt = ps.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK(lt > 0.078);
		CHECK(lt < 0.079);
	}

	SUBCASE("Coherent, late intersect.")
	{
		ps.is_clockwise = 1;
		ps.is_coherent = 1;

		auto particle = CreateParticle(0, 0, v2(7e-7, 4.3e-7), &ps);
		double lt = ps.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK(lt > 3.82); 
		CHECK(lt < 3.83);
	}

	SUBCASE("Repositioned.")
	{
		ps50.is_clockwise = 0;
		ps50.is_coherent = 1;

		auto particle = CreateParticle(0, 0, v2(6.5e-7, 1.35e-6), &ps50);
		double lt = ps50.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK(lt == 0); // Op dit moment betekend 0 dat een particle nergens op is gebotst.
	}

	SUBCASE("Repositioned.")
	{
		ps50.is_clockwise = 0;
		ps50.is_coherent = 1;

		auto particle = CreateParticle(0, 0, v2(6.15e-7, 1.35e-6), &ps50);
		double lt = ps50.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);
		
		CHECK(lt > 2.21);
		CHECK(lt < 2.215);
	}

	SUBCASE("Edge case 1.")
	{
		ps.is_clockwise = 0;
		ps.is_coherent  = 1;

		auto particle = CreateParticle(0, 1, v2(1.00000000e-06, 5.0000000e-07), &ps8);
		double lt = ps8.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK(lt == 0); // Op dit moment betekend 0 dat een particle nergens op is gebotst.
	}

	SUBCASE("Edge case 2.")
	{
		ps8.is_clockwise = 0;
		ps8.is_coherent  = 1;

		auto particle = CreateParticle(0, 1, v2(1.00000000e-06, 5.0000100e-07), &ps8);
		double lt = ps8.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK(lt < 1);
	}

	SUBCASE("Edge case 3.")
	{
		ps8.is_clockwise = 0;
		ps8.is_coherent  = 1;

		auto particle = CreateParticle(0, 1, v2(1.00000000e-06, 4.9999900e-07), &ps8);
		double lt = ps8.angular_speed * TraceOrbit(&particle, &impurity_settings, grid.GetImpurities(), grid.GetIndex(), pm);

		CHECK(lt < 1);
	}
}

