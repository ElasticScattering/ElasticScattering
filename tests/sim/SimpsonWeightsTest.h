#pragma once

#include "doctest.h"
#include "src/sim/es/util.h"

TEST_CASE("Weight1D test") {
	int dim = 7;

	double w = SimpsonWeight(0, dim);
	CHECK(w == 1);

	w = SimpsonWeight(dim - 1, dim);
	CHECK(w == 1);

	w = SimpsonWeight(1, dim);
	CHECK(w == 4);

	w = SimpsonWeight(dim - 2, dim);
	CHECK(w == 4);

	w = SimpsonWeight(dim - 3, dim);
	CHECK(w == 2);
}

TEST_CASE("Weight2D test") {
	int dim = 7;

	double w = SimpsonWeight2D(dim - 1, 0, dim);
	CHECK(w == 1);

	w = SimpsonWeight2D(0, 0, dim);
	CHECK(w == 1);

	w = SimpsonWeight2D(1, 0, dim);
	CHECK(w == 4);

	w = SimpsonWeight2D(2, 0, dim);
	CHECK(w == 2);

	w = SimpsonWeight2D(3, 0, dim);
	CHECK(w == 4);

	w = SimpsonWeight2D(0, 1, dim);
	CHECK(w == 4);

	w = SimpsonWeight2D(1, 1, dim);
	CHECK(w == 16);

	w = SimpsonWeight2D(3, 3, dim);
	CHECK(w == 16);

	w = SimpsonWeight2D(1, 2, dim);
	CHECK(w == 8);

	w = SimpsonWeight2D(2, 1, dim);
	CHECK(w == 8);

	w = SimpsonWeight2D(0, dim - 1, dim);
	CHECK(w == 1);
}
