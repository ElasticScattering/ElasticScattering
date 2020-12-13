#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/escl/cell_grid.h"

#include <vector>

void TestGridMovement(Orbit& o, Intersection& exit, std::vector<Intersection>& expected_intersections)
{
	Intersection entry;

	int cells_per_row = 10;
	v2 spawn_range = { -1e-6, 2e-6 };
	double cell_size = (spawn_range.y - spawn_range.x) / (double)cells_per_row;

	int i = 0;
	while (1)
	{
		entry = exit;
		bool cell_available = GetNextCell(&o, 0, cell_size, cells_per_row, spawn_range, &entry, &exit);

		if (!cell_available)
		{
			CHECK(i == (expected_intersections.size()));
			break;
		}

		CHECK_ALMOST(exit.position.x, expected_intersections[i].position.x);
		CHECK_ALMOST(exit.position.y, expected_intersections[i].position.y);
		CHECK_ALMOST(exit.incident_angle, expected_intersections[i].incident_angle);
		CHECK(exit.entering_cell.x == expected_intersections[i].entering_cell.x);
		CHECK(exit.entering_cell.y == expected_intersections[i].entering_cell.y);

		i++;
	}
}

TEST_CASE("GetNextCell")
{
	SUBCASE("1st setup.")
	{
		Orbit o({ 7e-07, 2.6858643180165056e-06 }, 2.3878643180165057e-06, true);

		Intersection exit;
		exit.position = v2(7e-7, 2.98e-7);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({ 5.000000000000001e-07, 3.063904263218665e-07}, {4,  4}, 6.199330215359630e+00),
			Intersection({ 2.000000000000000e-07, 3.509347684890618e-07}, {3,  4}, 6.072232037914050e+00),
			Intersection({-9.999999999999989e-08, 4.359985439632108e-07}, {2,  4}, 5.941550985064070e+00),
			Intersection({-2.611936248637251e-07, 5.000000000000000e-07}, {2,  5}, 5.868903316607379e+00),
			Intersection({-4.000000000000000e-07, 5.664549660374142e-07}, {1,  5}, 5.804443623619448e+00),
			Intersection({-7.000000000000000e-07, 7.514662200913574e-07}, {0,  5}, 5.656703934600270e+00),
			Intersection({-7.647224226072912e-07, 8.000000000000001e-07}, {0,  6}, 5.622823418415306e+00)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("2nd setup (same as first but with counterclockwise orbit.")
	{
		Orbit o({ 7e-07, -2.0898643180165058e-06 }, 2.3878643180165057e-06, false);

		Intersection exit;
		exit.position = v2(7e-7, 2.98e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({ 5.000000000000001e-07,  2.896095736781334e-07}, {4,  4}, 8.385509181995565e-02),
			Intersection({ 2.000000000000000e-07,  2.450652315109381e-07}, {3,  4}, 2.109532692655360e-01),
			Intersection({ 2.293471043684734e-08,  2.000000000000000e-07}, {3,  3}, 2.874880749540860e-01),
			Intersection({-9.999999999999989e-08,  1.600014560367891e-07}, {2,  3}, 3.416343221155165e-01),
			Intersection({-4.000000000000000e-07,  2.954503396258569e-08}, {1,  3}, 4.787416835601386e-01),
			Intersection({-6.199757562702194e-07, -9.999999999999989e-08}, {1,  2}, 5.857026755566084e-01),
			Intersection({-7.000000000000000e-07, -1.554662200913575e-07}, {0,  2}, 6.264813725793159e-01),
			Intersection({-9.870846416084292e-07, -4.000000000000000e-07}, {0,  1}, 7.845750310926594e-01)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("3rd setup.")
	{
		Orbit o({ 7e-07, 1.4919321590082528e-06 }, 1.1939321590082528e-06, true);

		Intersection exit;
		exit.position = v2(7.000000000000000e-07, 2.980000000000000e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({ 5.000000000000001e-07, 3.148705632053364e-07}, {4, 4}, 6.114878111366652e+00),
			Intersection({ 2.000000000000000e-07, 4.077393715787716e-07}, {3, 4}, 5.851079191355701e+00),
			Intersection({ 3.551178171517997e-08, 5.000000000000000e-07}, {3, 5}, 5.692952546289328e+00),
			Intersection({-9.999999999999989e-08, 6.056624481898653e-07}, {2, 5}, 5.548902669427402e+00),
			Intersection({-2.729870953123098e-07, 8.000000000000001e-07}, {2, 6}, 5.330553837030871e+00),
			Intersection({-4.000000000000000e-07, 1.027740387954405e-06}, {1, 6}, 5.111709509159294e+00),
			Intersection({-4.277691177937253e-07, 1.100000000000000e-06}, {1, 7}, 5.046860520666858e+00),
			Intersection({-4.903875328875842e-07, 1.400000000000000e-06}, {1, 8}, 4.789464754429844e+00),
			Intersection({-4.756622703221965e-07, 1.700000000000000e-06}, {1, 9}, 4.537223513199825e+00),
			Intersection({-4.000000000000000e-07, 1.956123930062101e-06}, {2, 9}, 4.313068451610086e+00)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("4th setup.")
	{
		Orbit o({ 7e-07 , -8.959321590082529e-07 }, 1.1939321590082528e-06, false);

		Intersection exit;
		exit.position = v2(7.000000000000000e-07, 2.980000000000000e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({ 5.000000000000001e-07,  2.811294367946634e-07}, {4,  4}, 1.683071958129345e-01),
			Intersection({ 2.262841535629006e-07,  2.000000000000000e-07}, {4,  3}, 4.079947636559396e-01),
			Intersection({ 2.000000000000000e-07,  1.882606284212283e-07}, {3,  3}, 4.321061158238857e-01),
			Intersection({-9.999999999999989e-08, -9.662448189865434e-09}, {2,  3}, 7.342826377521838e-01),
			Intersection({-1.899247151139072e-07, -9.999999999999989e-08}, {2,  2}, 8.410940816173249e-01),
			Intersection({-3.860595259816659e-07, -4.000000000000000e-07}, {2,  1}, 1.142438932907954e+00),
			Intersection({-4.000000000000000e-07, -4.317403879544051e-07}, {1,  1}, 1.171475798020293e+00),
			Intersection({-4.777455537510948e-07, -7.000000000000000e-07}, {1,  0}, 1.405944052499597e+00)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("5th setup.")
	{
		Orbit o({ 7e-07, 7.755728636033011e-07 }, 4.775728636033011e-07, true);

		Intersection exit;
		exit.position = v2(7.000000000000000e-07, 2.980000000000000e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({5.000000000000001e-07, 3.418957486315086e-07}, {4, 4}, 5.851079191355701e+00),
			Intersection({3.099545707283141e-07, 5.000000000000000e-07}, {4, 5}, 5.327473916493865e+00),
			Intersection({2.230522512290111e-07, 8.000000000000001e-07}, {4, 6}, 4.661218150259741e+00),
			Intersection({3.495390560708725e-07, 1.100000000000000e-06}, {4, 7}, 3.965546700212334e+00),
			Intersection({5.000000000000000e-07, 1.209249978575094e-06}, {5, 7}, 3.573698769413678e+00),
			Intersection({8.000000000000001e-07, 1.242558773508790e-06}, {6, 7}, 2.930639384324257e+00),
			Intersection({1.050460943929127e-06, 1.100000000000000e-06}, {6, 6}, 2.317638606967253e+00),
			Intersection({1.100000000000000e-06, 1.036486335946638e-06}, {7, 6}, 2.148775272056415e+00),
			Intersection({1.176947748770989e-06, 8.000000000000002e-07}, {7, 5}, 1.621967156919846e+00),
			Intersection({1.100000000000000e-06, 5.146593912599645e-07}, {6, 5}, 9.928173815333787e-01),
			Intersection({1.090045429271686e-06, 5.000000000000001e-07}, {6, 4}, 9.557113906857211e-01),
			Intersection({8.000000000000002e-07, 3.085869536978123e-07}, {5, 4}, 2.109532692655360e-01)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("6th setup.")
	{
		Orbit o({ 7e-07,  -1.7957286360330114e-07 }, 4.775728636033011e-07, false);

		Intersection exit;
		exit.position = v2(7.000000000000000e-07, 2.980000000000000e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({5.000000000000001e-07,  2.541042513684914e-07}, {4, 4}, 4.321061158238857e-01),
			Intersection({4.101719798462421e-07,  2.000000000000000e-07}, {4, 3}, 6.521254451677465e-01),
			Intersection({2.291029842648950e-07, -9.999999999999989e-08}, {4, 2}, 1.403396270295512e+00),
			Intersection({2.763400920665157e-07, -4.000000000000000e-07}, {4, 1}, 2.050545952078872e+00),
			Intersection({5.000000000000000e-07, -6.132499785750936e-07}, {5, 1}, 2.709486537765907e+00),
			Intersection({8.000000000000001e-07, -6.465587735087899e-07}, {6, 1}, 3.352545922855330e+00),
			Intersection({1.100000000000000e-06, -4.404863359466380e-07}, {7, 1}, 4.134410035123171e+00),
			Intersection({1.123659907933484e-06, -4.000000000000000e-07}, {7, 2}, 4.232639355100714e+00),
			Intersection({1.170897015735105e-06, -1.000000000000000e-07}, {7, 3}, 4.879789036884074e+00),
			Intersection({1.100000000000000e-06,  8.134060874003542e-08}, {6, 3}, 5.290367925646208e+00),
			Intersection({9.898280201537577e-07,  2.000000000000001e-07}, {6, 4}, 5.631059862011840e+00),
			Intersection({8.000000000000002e-07,  2.874130463021876e-07}, {5, 4}, 6.072232037914050e+00)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("7th setup.")
	{
		Orbit o({ 7e-07, 4.1170782466745266e-07 }, 1.1370782466745266e-07, true);

		Intersection exit;
		exit.position = v2(7.000000000000000e-07, 2.980000000000000e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({6.283483345228397e-07, 5.000000000000000e-07}, {5,  5}, 3.823324077988501e+00),
			Intersection({7.716516654771601e-07, 5.000000000000001e-07}, {5,  4}, 2.459861229191086e+00),
			Intersection({8.000000000000001e-07, 4.658323975581844e-07}, {6,  4}, 2.066893572964171e+00),
			Intersection({8.000000000000002e-07, 3.575832517767211e-07}, {5,  4}, 1.074699080625624e+00)
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("8th setup.")
	{
		Orbit o({ 7e-07, 3.9351457272066023e-07 }, 9.551457272066024e-08, true);

		Intersection exit;
		exit.position = v2();
		exit.dphi = PI2;
		exit.entering_cell = { 5, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
		};

		std::vector<v2i> expected_cells{
		};

		std::vector<double> expected_cell_phis{
		};

		TestGridMovement(o, exit, intersections);
	}

	SUBCASE("9th setup.")
	{
		Orbit o({ 4.623049896354275e-07, -1.3769501036457238e-07 }, 4.775728636033011e-07, false);

		Intersection exit;
		exit.position = v2(4.623049896354275e-07, 3.398778532387288e-07);
		exit.dphi = PI2;
		exit.entering_cell = { 4, 4 };
		exit.incident_angle = 0;

		std::vector<Intersection> intersections{
			Intersection({ 2.000000000000000e-07,  2.613938675343432e-07}, {3, 4}, 5.814617014956047e-01),
			Intersection({ 1.246099792708551e-07,  2.000000000000000e-07}, {3, 3}, 7.853981633974483e-01),
			Intersection({-1.377790886426873e-08, -9.999999999999989e-08}, {3, 2}, 1.491783754263309e+00),
			Intersection({ 6.321611173651196e-08, -4.000000000000000e-07}, {3, 1}, 2.152258028290502e+00),
			Intersection({ 2.000000000000001e-07, -5.367838882634881e-07}, {4, 1}, 2.560130952094189e+00),
			Intersection({ 5.000000000000000e-07, -6.137779088642687e-07}, {5, 1}, 3.220605226121379e+00),
			Intersection({ 8.000000000000001e-07, -4.753900207291446e-07}, {6, 1}, 3.926990816987241e+00),
			Intersection({ 8.613938675343431e-07, -4.000000000000000e-07}, {6, 2}, 4.130927278889084e+00),
			Intersection({ 9.383878881351238e-07, -1.000000000000000e-07}, {6, 3}, 4.791401552916278e+00),
			Intersection({ 8.000000000000002e-07,  1.999999999999997e-07}, {5, 4}, 5.497787143782137e+00),
			Intersection({ 5.000000000000001e-07,  3.383878881351238e-07}, {4, 4}, 6.204172734647999e+00)
		};

		TestGridMovement(o, exit, intersections);
	}
}

/* Todo: check alle waardes van de intersectie? */
/*
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
		CHECK_APPROX(expected_intersection.x, i.position.x);
		CHECK_APPROX(expected_intersection.y, i.position.y);
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
*/