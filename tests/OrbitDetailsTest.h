#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/sim/es/details.h"

TEST_CASE("InsideImpurity")
{
	SUBCASE("returns true")
	{
		bool inside = InsideImpurity({ 1, 1 }, { 0, 0 }, 2);
		CHECK(inside == true);
	}

	SUBCASE("returns true")
	{
		bool inside = InsideImpurity({ 0, 0 }, { 1, 1 }, 2);
		CHECK(inside == true);
	}

	SUBCASE("returns false")
	{
		bool inside = InsideImpurity({ 2, 1 }, { 0, 0 }, 2);
		CHECK(inside == false);
	}

	SUBCASE("Position on edge of impurity returns false")
	{
		bool inside = InsideImpurity({ 1, 0 }, { 0, 0 }, 1);
		CHECK(inside == false);
	}
}

TEST_CASE("Cyclotron Orbit")
{
	v2 pos = { 1e-6, 1e-6 };
	v2 vel = { 7e5, 0 };
	double wc = E * 2.5 / M0;
	double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
	double radius = vf / wc;

	SUBCASE("Electron is shifted up in y with positive vx")
	{
		auto c = GetCyclotronOrbitCenter(pos, vel, radius, vf, true);
		CHECK(c.x == pos.x);
		CHECK(c.y == pos.y + radius);
	}

	SUBCASE("Hole is shifted down in y with positive vx")
	{
		auto c = GetCyclotronOrbitCenter(pos, vel, radius, vf, false);
		CHECK(c.x == pos.x);
		CHECK(c.y == pos.y - radius);
	}
}

TEST_CASE("BoundTime Direction")
{
	double phi = -0.25;
	double alpha = 1;
	double w = 0.5;

	SUBCASE("Clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, false, true);
		CHECK(t == 0.25);
	}

	SUBCASE("Counter-is_clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, true, true);
		CHECK(t == 0.75);
	}

	SUBCASE("Clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, true, false);
		CHECK(t == 0.25);
	}

	SUBCASE("Counter-is_clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, false, false);
		CHECK(t == 0.75);
	}
}

TEST_CASE("BoundTime Sectors")
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(-0.25, w, alpha, false, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(PI / 2 - 0.25, w, alpha, false, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(PI - 0.25, w, alpha, false, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(3 * PI / 2 - 0.25, w, alpha, false, false, true);
	CHECK(t == 0.25);
}

TEST_CASE("BoundTime Corners")
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(0.5, w, alpha, false, false, true);
	CHECK(t == 1);

	t = GetBoundTime(0.5, w, alpha, false, false, false);
	CHECK(t == 0);
}

TEST_CASE("Circles Cross")
{
	Orbit orbit1({ 1, 2 }, 1, true, 0, 0);
	CHECK(!CirclesCross(&orbit1, { 5, 5 }, 1));
	CHECK(!CirclesCross(&orbit1, { 1, 3 }, 5));
	
	Orbit orbit2({ 0, 0 }, 1, true, 0, 0);
	CHECK(CirclesCross(&orbit2, { 0, 0.5 }, 0.6));
	
	Orbit orbit3({ 0, 0 }, 1000, true, 0, 0);
	CHECK(CirclesCross(&orbit3, { 2, 0 }, 1001));
}

TEST_CASE("Circle Crosspoints")
{
	SUBCASE("Symmetric")
	{
		Orbit orbit1({ -1, 0 }, 1.5, true, 0, 0);

		auto ps = GetCrossPoints(&orbit1, { 1, 0 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };
		CHECK(p1.x == 0);
		CHECK(abs(p1.y) == sqrt(1.5 * 1.5 - 1));

		CHECK(p2.x == 0);
		CHECK(abs(p2.y) == sqrt(1.5 * 1.5 - 1));
	}

	SUBCASE("Somewhat Symmetric")
	{
		Orbit orbit1({ 0, 0 }, 1, true, 0, 0);

		auto ps = GetCrossPoints(&orbit1, { 1, 1 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };
		CHECK_APPROX(1, pow(p1.x, 2) + pow(p1.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x - 1, 2) + pow(p1.y - 1, 2));

		CHECK_APPROX(1, pow(p2.x, 2) + pow(p2.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x - 1, 2) + pow(p2.y - 1, 2));
	}

	SUBCASE("Asymmetric")
	{
		Orbit orbit1({ 100, 0 }, 100, true, 0, 0);

		auto ps = GetCrossPoints(&orbit1, { -1, 0 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };

		CHECK_APPROX(pow(100, 2), pow(p1.x - 100, 2) + pow(p1.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x + 1, 2) + pow(p1.y, 2));

		CHECK_APPROX(pow(100, 2), pow(p2.x - 100, 2) + pow(p2.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x + 1, 2) + pow(p2.y, 2));
	}

	SUBCASE("Asymmetric, inside")
	{
		Orbit orbit1({ 100, -0.5 }, 100, true, 0, 0);

		auto ps = GetCrossPoints(&orbit1, { 1, 0.5 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };

		CHECK_APPROX(pow(100, 2), pow(p1.x - 100, 2) + pow(p1.y + 0.5, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x - 1, 2) + pow(p1.y - 0.5, 2));

		CHECK_APPROX(pow(100, 2), pow(p2.x - 100, 2) + pow(p2.y + 0.5, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x - 1, 2) + pow(p2.y - 0.5, 2));
	}
}

TEST_CASE("Cross Time")
{
	const v2 pos = { 3, 4 };

	Orbit orbit1({ 0,0 }, 5, true, 0, 0);
	Orbit orbit2({ 0,0 }, 5, false, 0, 0);

	orbit1.particle_angle = GetAngle(pos, &orbit1);
	orbit2.particle_angle = GetAngle(pos, &orbit2);

	double ir = 0.1;
	double w = 2;

	double t = GetFirstCrossTime(&orbit1, { 5, 0 }, ir, w, { 0, 0});
	CHECK_APPROX_LOW(t, 0.907 / 2);

	double t2 = GetFirstCrossTime(&orbit2, { 5, 0 }, ir, w, { -0, 0 });
	CHECK_APPROX(t + t2, (PI2 - (ir * 2.0 / orbit2.radius)) / w);

	ir = 0.059;
	w = 100;

	t = GetFirstCrossTime(&orbit1, { 5, 0 }, ir, w, { 0, 0 });
	t2 = GetFirstCrossTime(&orbit2, { 5, 0 }, ir, w, { 0, 0 });
	CHECK_APPROX(t + t2, (PI2 - (ir * 2.0 / orbit1.radius)) / w);
}
