#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/escl/cell_grid.h"

TEST_CASE("GetNextCell")
{
	Orbit o({ 1.5, 1.5 }, sqrt(2.5), true);
	Intersection entry;
	entry.position = v2(1, 0);
	entry.dphi = GetAngle({ 1, 0 }, &o);
	entry.entering_cell = { 0, 0 };

	Intersection exit;

	//xlow, ylow, L, xc, yc, rc, isEl, entry):
	bool cell_available = GetNextCell(&o, entry, 3, 10, { 0, 30 }, &exit);

}

TEST_CASE("UpdateBestIntersect Tests")
{
	Orbit o({ 1.5, 1.5 }, sqrt(2.5), true);
	Intersection entry;
	entry.position = v2(1, 0);
	entry.dphi = GetAngle({ 1, 0 }, &o);

	Intersection candidate;
	Intersection closest;

	UpdateBestIntersect(candidate, { 0,1 }, entry, true, 20, 1e-2, &closest);
}

/* Todo: check alle waardes van de intersectie? */
TEST_CASE("GetFirstBoundaryIntersects Tests")
{
	SUBCASE("Hit 1")
	{
		Orbit orbit({ 0, 0 }, 2, false, 0, 0);
		Intersection i;

		bool hit = GetFirstBoundaryIntersect(&orbit, { 0, -3 }, { 0, 3 }, 6, 0, &i);

		v2 expected_intersection = { 0, -2 };
		CHECK(hit);
		CHECK(expected_intersection == i.position);
	}

	SUBCASE("Hit 2")
	{
		Orbit orbit({ 0.7, 0 }, .5, false, 0, 0);
		Intersection i;

		bool hit = GetFirstBoundaryIntersect(&orbit, { 1, 0 }, { 1, 6 }, 6, 0, &i);

		v2 expected_intersection = { 1, 0.4 };
		CHECK(hit);
		CHECK_APPROX(expected_intersection.x, i.position.x, "");
		CHECK_APPROX(expected_intersection.y, i.position.y, "");
	}

	SUBCASE("Hit 3")
	{
		Orbit orbit({ 1, 1 }, 1.1, false, 0, 0);
		Intersection i;

		bool hit = GetFirstBoundaryIntersect(&orbit, { 0, 0 }, { 0, 2 }, 2, 0, &i);

		v2 expected_intersection = { 0, (1 - sqrt(1.1 * 1.1 - 1)) };
		CHECK(hit);
		CHECK((expected_intersection == i.position));
	}

	SUBCASE("Orbit that encircles line segment should return no intersection")
	{
		Orbit orbit({ 0.2 , 0.2 }, 20, false, 0, 0);

		Intersection i;
		bool hit = GetFirstBoundaryIntersect(&orbit, { 0, -3 }, { 0, 3 }, 6, 0, &i);

		CHECK(!hit);
	}

	SUBCASE("Orbit that misses the line segment should return no intersection")
	{
		Orbit orbit({ 0, 0 }, 1, false, 0, 0);

		Intersection i;
		bool hit = GetFirstBoundaryIntersect(&orbit, { 5, -50 }, { 5, 60 }, 110, 0, &i);

		CHECK(!hit);
	}
}


TEST_CASE("get_cell Tests")
{
	SUBCASE("Extends is empty, (0,0) should get the first cell") {
		auto cell = get_cell(0, 0, { 0, 0.4 }, 6);
		
		auto expected_cell = v2i(0, 0);
		CHECK(cell == expected_cell);
	}

	SUBCASE("Extends is not empty, (0,0) should be offset") {
		auto cell = get_cell(0, 0, { -0.1, 0.4 }, 6);
		
		auto expected_cell = v2i(1, 1);
		CHECK(cell == expected_cell);
	}

	SUBCASE("Position in the middle should be in the middle cell.") {
		int cells_per_row = 9;
		v2 range = { -0.1, 0.4 };
		double dim = range.y - range.x;
		double mid = range.x + 0.5 * dim;
		auto cell = get_cell(mid, mid, range, cells_per_row);

		v2i expected_cell = { cells_per_row / 2, cells_per_row / 2 };
		CHECK(cell == expected_cell);
	}

	SUBCASE("Position at the bottom right corner should be in the last cell.") {
		int cells_per_row = 9;
		auto cell = get_cell(0.399, 0.399, { -0.1, 0.4 }, cells_per_row);
		
		v2i expected_cell = { cells_per_row-1, cells_per_row-1 };
		CHECK(cell == expected_cell);
	}
}


TEST_CASE("get_cell_index Tests")
{
	SUBCASE("Extends is empty, (0,0) should get the first cell.") {
		auto cell = get_cell_index({ 0, 0 }, { 0, 0.4 }, 6);
		CHECK(cell == 0);
	}

	SUBCASE("Extends is not empty, (0,0) should not be the first cell.") {
		auto cell = get_cell_index({ 0, 0 }, { -0.1, 0.4 }, 6);
		CHECK(cell == 7);
	}

	SUBCASE("Position at the bottom right corner should be in the last cell.") {
		int cells_per_row = 9;
		auto cell = get_cell_index({ 0.399, 0.399 }, { -0.1, 0.4 }, cells_per_row);
		
		int last_cell = (cells_per_row * cells_per_row) - 1;
		CHECK(cell == last_cell);
	}
}


TEST_CASE("to_index Tests")
{
	CHECK(get_index({ 0, 0 }, 30) ==  0);
	CHECK(get_index({ 3, 2 }, 30) == 63);
	CHECK(get_index({ 2, 3 }, 30) == 92);
	CHECK(get_index({ 4, 3 }, 10) == 34);
}


TEST_CASE("within_bounds Tests")
{
	SUBCASE("Negative position X should be false.")
	{
		bool result = within_bounds({ -1, 3 }, 30);
		CHECK(result == false);
	}

	SUBCASE("Negative position Y should be false.")
	{
		bool result = within_bounds({ 3, -3 }, 30);
		CHECK(result == false);
	}

	SUBCASE("Out of bounds position X should be false.")
	{
		bool result = within_bounds({ 30, 3 }, 30);
		CHECK(result == false);
	}

	SUBCASE("Out of bounds position Y should be false.")
	{
		bool result = within_bounds({ 30, 12 }, 30);
		CHECK(result == false);
	}

	SUBCASE("First cell should be true.")
	{
		bool result = within_bounds({ 0, 0 }, 30);
		CHECK(result == true);
	}

	SUBCASE("Cell in range should be true.")
	{
		bool result = within_bounds({ 12, 12 }, 30);
		CHECK(result == true);
	}

	SUBCASE("Last cell should be true.")
	{
		bool result = within_bounds({ 29, 29 }, 30);
		CHECK(result == true);
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
