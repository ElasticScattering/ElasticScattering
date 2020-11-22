#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/escl/impurity_grid.h"

TEST_CASE("GetFirstBoundaryIntersects Tests")
{
	SUBCASE("1")
	{
		Orbit orbit({ 0,0 }, 2, false, 0, 0);
		double L = 6;
		Intersection i;
		v2 expected_intersection = { 0, 2 };

		bool hit = GetFirstBoundaryIntersect({ 0, -3 }, { 0, 3 }, L, &orbit, 0, &i);

		CHECK(hit == true);
		CHECK((expected_intersection == i.position));
	}

	SUBCASE("2")
	{
		Orbit orbit({ 0.7, 0 }, .5, false, 0, 0);
		double L = 6;
		Intersection i;
		v2 expected_intersection = { 1, 0.4 };

		bool hit = GetFirstBoundaryIntersect({ 1, 0 }, { 1, 6 }, L, &orbit, 0, &i);

		CHECK(hit == true);
		CHECK((expected_intersection == i.position));
	}

	SUBCASE("3")
	{
		Orbit orbit({ 1, 1 }, 1.1, false, 0, 0);
		double L = 2;
		Intersection i;
		v2 expected_intersection = { 0, (1 - sqrt(1.1*1.1 -1)) };

		bool hit = GetFirstBoundaryIntersect({ 0, 0 }, { 0, 2 }, L, &orbit, 0, &i);

		CHECK(hit == true);
		CHECK((expected_intersection == i.position));
	}

	SUBCASE("Orbit that encircles line segment should return no intersection")
	{
		Orbit orbit({ 0.2 , 0.2 }, 20, false, 0, 0);
		double L = 6;

		Intersection i;
		bool hit = GetFirstBoundaryIntersect({ 0, -3 }, { 0, 3 }, L, &orbit, 0, &i);

		CHECK(hit == false);
	}

	SUBCASE("Orbit that misses the line segment should return no intersection")
	{
		Orbit orbit({ 0, 0 }, 1, false, 0, 0);
		double L = 110;

		Intersection i;
		bool hit = GetFirstBoundaryIntersect({ 5, -50 }, { 5, 60 }, L, &orbit, 0, &i);

		CHECK(hit == false);
	}
}

TEST_CASE("to_grid Tests")
{
	SUBCASE("Extends is empty, (0,0) should get the first cell") {
		auto cell = to_grid(0, 0, { 0, 0.4 }, 6);
		auto expected_cell = v2i(0, 0);
		CHECK(cell == expected_cell);
	}

	SUBCASE("Extends is not empty, (0,0) should be offset") {
		auto cell = to_grid(0, 0, { -0.1, 0.4 }, 6);
		auto expected_cell = v2i(1, 1);
		CHECK(cell == expected_cell);
	}

}


/*
TEST_CASE("GetNextCell Tests")
{
	Orbit orbit({ 0, 0 }, 1, false, 0, 0);

	Intersection start_intersect;
	start_intersect.position = { cos(0), sin(0) };
	start_intersect.dphi = 0;

	Intersection i;

	// GetNextCell is wat groter dan de corresponderende python test (find_exit_intersect)
	// Iig setup zorgen die tot dezelfde celzijde intersecties leiden.
	//auto result = GetNextCell(orbit, , start_intersect, 0.8, 30, &i);
}
*/

//boxmotion, weet niet of dit hoeft voor mij.
