#pragma once

#include <doctest.h>
#include "src/scattering/escl/cell_grid.h"

TEST_CASE("DifferentPoint")
{
	v2 p1 = { 5.000000e-07, 3.063904e-07 };

	CHECK(!DifferentPoint(p1, p1, 3e-7));
	CHECK(!DifferentPoint(p1, { 5.000000e-07, 3.063904e-07 }, 3e-7));

	v2 p3 = { 2.000000e-07, 3.063904e-07 };
	CHECK(DifferentPoint(p1, p3, 3e-7));
}


TEST_CASE("to_world Tests")
{
	SUBCASE("Setup from GetNextCell test.")
	{
		int cells_per_row = 10;
		v2 spawn_range = { -1e-6, 2e-6 };
		double start = spawn_range.x;
		double size = spawn_range.y - spawn_range.x;
		double cell_size = size / (double)cells_per_row;

		v2 pos = to_world({ 5, 4 }, start, cell_size);

		v2 expected_pos = {
			start + 5 / (double)cells_per_row * size,
			start + 4 / (double)cells_per_row * size
		};

		CHECK(pos.x == expected_pos.x);
		CHECK_ALMOST(pos.y, expected_pos.y);
	}

	SUBCASE("Rounded numbers.")
	{
		int cells_per_row = 10;
		v2 spawn_range = { -10, 20 };
		double start = spawn_range.x;
		double size = spawn_range.y - spawn_range.x;
		double cell_size = size / (double)cells_per_row;

		v2 pos = to_world({ 2, 3 }, start, cell_size);

		v2 expected_pos = v2(spawn_range.x + 2 * cell_size, spawn_range.x + 3 * cell_size);

		CHECK(pos.x == expected_pos.x);
		CHECK(pos.y == expected_pos.y);
	}

	SUBCASE("Rounded numbers at origin.")
	{
		int cells_per_row = 10;
		v2 spawn_range = { -10, 20 };
		double start = spawn_range.x;
		double size = spawn_range.y - spawn_range.x;
		double cell_size = size / (double)cells_per_row;

		v2 pos = to_world({ 0, 0 }, start, cell_size);

		v2 expected_pos = { spawn_range.x, spawn_range.x };

		CHECK(pos.x == expected_pos.x);
		CHECK(pos.y == expected_pos.y);
	}
}


TEST_CASE("get_cell Tests")
{
	SUBCASE("Extends is empty, (0,0) should get the first cell") {
		auto cell = get_cell(v2(0, 0), 0, 0.4, 6);

		auto expected_cell = v2i(0, 0);
		CHECK(cell == expected_cell);
	}

	SUBCASE("Extends is not empty, (0,0) should be offset") {
		auto cell = get_cell(v2(0, 0), -0.1, 0.4 + 0.1, 6);

		auto expected_cell = v2i(1, 1);
		CHECK(cell == expected_cell);
	}

	SUBCASE("Position in the middle should be in the middle cell.") {
		int cells_per_row = 9;
		v2 range = { -0.1, 0.4 };
		double dim = range.y - range.x;
		double mid = range.x + 0.5 * dim;
		auto cell = get_cell(v2(mid, mid), range.x, dim, cells_per_row);

		v2i expected_cell = { cells_per_row / 2, cells_per_row / 2 };
		CHECK(cell == expected_cell);
	}

	SUBCASE("Position at the bottom right corner should be in the last cell.") {
		int cells_per_row = 9;
		auto cell = get_cell(v2(0.399, 0.399), -0.1, 0.4 + 0.1, cells_per_row);

		v2i expected_cell = { cells_per_row - 1, cells_per_row - 1 };
		CHECK(cell == expected_cell);
	}
}

TEST_CASE("to_index Tests")
{
	CHECK(get_index({ 0, 0 }, 30) == 0);
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
