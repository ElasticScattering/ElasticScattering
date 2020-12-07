#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/escl/details.h"
#include <stdio.h>


TEST_CASE("GetAngle / AngleVelocity")
{
	double angle;

	SUBCASE("Clockwise")
	{
		Orbit o({ 0, 0 }, 1, true);

		angle = GetAngle({ 1, 0 }, &o);
		CHECK_APPROX(angle, HALF_PI);

		angle = GetAngle({ -1, 0 }, &o);
		CHECK_APPROX(angle, (3 * HALF_PI));

		angle = GetAngle({ 0, 1 }, &o);
		CHECK_APPROX(angle, PI);

		angle = GetAngle({ 0, -1 }, &o);
		CHECK_APPROX(angle, 0);
	}

	SUBCASE("Counterclockwise")
	{
		Orbit o({ 0, 0 }, 1, false);

		angle = GetAngle({ -1, 0 }, &o);
		CHECK_APPROX(angle, HALF_PI);

		angle = GetAngle({ 0, 1 }, &o);
		CHECK_APPROX(angle, 0);

		angle = GetAngle({ 1, 0 }, &o);
		CHECK_APPROX(angle, (3 * HALF_PI));

		angle = GetAngle({ 0, -1 }, &o);
		CHECK_APPROX(angle, PI);
	}

	SUBCASE("Clockwise #2")
	{
		Orbit o({ 6, 6 }, sqrt(2), true);
		angle = GetAngle({ 5, 5 }, &o);
		CHECK_APPROX(angle, (7 * PI / 4));
	}

	SUBCASE("Counterclockwise #2")
	{
		Orbit o({ 6, 6 }, sqrt(2), false);
		angle = GetAngle({ 5, 5 }, &o);
		CHECK_APPROX(angle, (3 * PI / 4));
	}
}

TEST_CASE("AngleInRange")
{
	SUBCASE("Direct inequality")
	{
		bool result = AngleInRange(1, { 0.5, 1.5 }, false);
		CHECK(result);
	}

	SUBCASE("")
	{
		bool result = AngleInRange(1, { 0.5, 1.5 }, true);
		CHECK(!result);
	}

	SUBCASE("")
	{
		bool result = AngleInRange(2, { 0.5, 1.5 }, true);
		CHECK(result);
	}

	SUBCASE("")
	{
		bool result = AngleInRange(0.1, { 0.5, 1.5 }, true);
		CHECK(result);
	}

	SUBCASE("Go up from 1.5 through 2pi to 0.5")
	{
		bool result = AngleInRange(1, { 1.5, 0.5 }, false);
		CHECK(!result);
	}

	SUBCASE("Go down from 1.to 0.5, direct inequality.")
	{
		bool result = AngleInRange(1, { 1.5, 0.5 }, true);
		CHECK(result);
	}

	SUBCASE("Full circle.")
	{
		bool result = AngleInRange(1, { 1, 1 }, true);
		CHECK(result);

		result = AngleInRange(2, { 1, 1 }, true);
		CHECK(result);
	}
}

TEST_CASE("Cross Angle")
{
	double a = GetCrossAngle(6.1, 0.1, true);
	CHECK(a == 6);

	a = GetCrossAngle(6.1, 0.1, false);
	CHECK_APPROX(a, PI2 - 6.0);
}
