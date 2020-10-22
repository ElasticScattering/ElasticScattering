#pragma once

#include "doctest.h"
#include "src/escl/weights.h"

TEST_CASE("Weight1D test") {
	int dim = 7;

	double w = GetWeight(0, dim);
	CHECK(w == 1);

	w = GetWeight(dim - 1, dim);
	CHECK(w == 1);

	w = GetWeight(1, dim);
	CHECK(w == 4);

	w = GetWeight(dim - 2, dim);
	CHECK(w == 4);

	w = GetWeight(dim - 3, dim);
	CHECK(w == 2);
}

TEST_CASE("Weight2D test") {
	int dim = 7;

	double w = GetWeight2D(dim - 1, 0, dim);
	CHECK(w == 1);

	w = GetWeight2D(0, 0, dim);
	CHECK(w == 1);

	w = GetWeight2D(1, 0, dim);
	CHECK(w == 4);

	w = GetWeight2D(2, 0, dim);
	CHECK(w == 2);

	w = GetWeight2D(3, 0, dim);
	CHECK(w == 4);

	w = GetWeight2D(0, 1, dim);
	CHECK(w == 4);

	w = GetWeight2D(1, 1, dim);
	CHECK(w == 16);

	w = GetWeight2D(3, 3, dim);
	CHECK(w == 16);

	w = GetWeight2D(1, 2, dim);
	CHECK(w == 8);

	w = GetWeight2D(2, 1, dim);
	CHECK(w == 8);

	w = GetWeight2D(0, dim - 1, dim);
	CHECK(w == 1);
}
