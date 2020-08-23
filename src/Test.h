#ifndef TEST_H
#define TEST_H

#include "ElasticScattering.h"
#include "Details.h"
#include <assert.h>
#include "doctest.h"

#define EPSILON 0.0001
#define EPSILON_HIGH 0.01
#define CHECK_APPROX(a, b)  { CHECK(abs((a)-(b)) < EPSILON); }
#define CHECK_APPROX_HIGH(a, b)  { assert(abs((a)-(b)) < EPSILON_HIGH); }

#define CHECK_CPU_GPU                                                            \
	e->Init(ip, sp);                                                              \
	double cpu_result = e->Compute();                                            \
	e2->Init(ip, sp);                                                             \
	double gpu_result = e2->Compute();                                           \
	std::cout << "CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << abs(gpu_result-cpu_result) << std::endl;  \
	CHECK_APPROX(cpu_result, gpu_result);  

TEST_CASE("Average lifetime on CPU and GPU")
{
	SimulationParameters sp;
	sp.region_size = 1e-6;
	sp.dim = 128;
	sp.particle_count = sp.dim * sp.dim;
	sp.particle_speed = 7e5;
	sp.particle_mass = 5 * M0;
	sp.impurity_count = 100;
	sp.impurity_radius = 1.5e-8;
	sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
	sp.alpha = PI / 4.0;
	sp.phi = 0;
	sp.magnetic_field = 0;
	sp.angular_speed = E * sp.magnetic_field / sp.particle_mass;
	sp.tau = 1e-12;

	InitParameters ip;
	ip.mode = Mode::AVG_LIFETIME;
	ip.show_info = false;
	ip.run_tests = true;

	auto e = new CPUElasticScattering();
	auto e2 = new GPUElasticScattering();

	SUBCASE("Standard test") {
		CHECK_CPU_GPU
	}

	SUBCASE("More impurities") {
		sp.particle_count = sp.dim * sp.dim;
		sp.impurity_count = 1000;
		
		CHECK_CPU_GPU;
	}
	
	SUBCASE("Larger impurities")
	{
		sp.impurity_radius = 1.5e-8;
		sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;

		CHECK_CPU_GPU
	}
	
	SUBCASE("With magnetic field")
	{
		sp.magnetic_field = 0;

		CHECK_CPU_GPU
	}

	SUBCASE("Different angle")
	{
		sp.phi = -sp.alpha - 1e-10;
		
		CHECK_CPU_GPU
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
		auto c = GetCyclotronOrbit(pos, vel, radius, vf, true);
		CHECK(c.x == pos.x);
		CHECK(c.y == pos.y + radius);
	}
	
	SUBCASE("Hole is shifted down in y with positive vx")
	{
		auto c = GetCyclotronOrbit(pos, vel, radius, vf, false);
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
		double t = GetBoundTime(phi, w, alpha, false, true);
		CHECK(t == 0.25);
	}

	SUBCASE("Counter-clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, true, true);
		CHECK(t == 0.75);
	}

	SUBCASE("Clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, true, false);
		CHECK(t == 0.25);
	}

	SUBCASE("Counter-clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, false);
		CHECK(t == 0.75);
	}
}

TEST_CASE("BoundTime Sectors")
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(-0.25, w, alpha, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(PI / 2 - 0.25, w, alpha, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(PI - 0.25, w, alpha, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(3 * PI / 2 - 0.25, w, alpha, false, true);
	CHECK(t == 0.25);
}

TEST_CASE("BoundTime Corners")
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(0.5, w, alpha, false, true);
	CHECK(t == 1);

	t = GetBoundTime(0.5, w, alpha, false, false);
	CHECK(t == 0);
}

TEST_CASE("Circles Cross")
{
	CHECK(!CirclesCross({ 1, 2 }, 1,    { 5, 5 }, 1));
	CHECK(!CirclesCross({ 1, 2 }, 1,    { 1, 3 }, 5));
	CHECK( CirclesCross({ 0, 0 }, 1,    { 0, 0.5 }, 0.6));
	CHECK( CirclesCross({ 0, 0 }, 1000, { 2, 0 }, 1001));
}

TEST_CASE("Circle Crosspoints")
{
	SUBCASE("Symmetric")
	{
		auto ps = GetCrossPoints({ -1, 0 }, 1.5, { 1, 0 }, 1.5);
		v2 p1 = ps.first;
		v2 p2 = ps.second;
		CHECK(p1.x == 0);
		CHECK(abs(p1.y) == sqrt(1.5 * 1.5 - 1));

		CHECK(p2.x == 0);
		CHECK(abs(p2.y) == sqrt(1.5 * 1.5 - 1));
	}
	
	SUBCASE("Somewhat Symmetric")
	{
		auto ps = GetCrossPoints({ 0, 0 }, 1, { 1, 1 }, 1.5);
		v2 p1 = ps.first;
		v2 p2 = ps.second;
		CHECK_APPROX(1,			  pow(p1.x, 2)     + pow(p1.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x - 1, 2) + pow(p1.y - 1, 2));

		CHECK_APPROX(1,			  pow(p2.x, 2)	   + pow(p2.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x - 1, 2) + pow(p2.y - 1, 2));
	}

	SUBCASE("Asymmetric")
	{
		auto ps = GetCrossPoints({ 100, 0 }, 100, { -1, 0 }, 1.5);
		v2 p1 = ps.first;
		v2 p2 = ps.second;

		CHECK_APPROX(pow(100, 2), pow(p1.x - 100, 2) + pow(p1.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x + 1, 2)   + pow(p1.y, 2));

		CHECK_APPROX(pow(100, 2), pow(p2.x - 100, 2) + pow(p2.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x + 1, 2)   + pow(p2.y, 2));
	}

	SUBCASE("Asymmetric, inside")
	{
		auto ps = GetCrossPoints({ 100, -0.5 }, 100, { 1, 0.5 }, 1.5);
		v2 p1 = ps.first;
		v2 p2 = ps.second;
		CHECK_APPROX(pow(100, 2), pow(p1.x - 100, 2) + pow(p1.y + 0.5, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x - 1, 2)   + pow(p1.y - 0.5, 2));

		CHECK_APPROX(pow(100, 2), pow(p2.x - 100, 2) + pow(p2.y + 0.5, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x - 1, 2)   + pow(p2.y - 0.5, 2));
	}
}

TEST_CASE("Phi")
{
	const v2 p1 = { 0, 0 };

	double phi = GetPhi({ 1, 0 }, p1, 1);
	CHECK(phi == 0);

	phi = GetPhi({ 0, 1 }, p1, 1);
	CHECK(phi == PI / 2);

	phi = GetPhi({ -1, 0 }, p1, 1);
	CHECK(phi == PI);

	phi = GetPhi({ 0, -1 }, p1, 1);
	CHECK(phi == 3 * PI / 2);

	phi = GetPhi({ 0.99999, -0.0045 }, p1, 1);
	CHECK_APPROX_HIGH(phi, PI2);
}

TEST_CASE("Cross Angle")
{
	double a = GetCrossAngle(6.1, 0.1, true);
	CHECK(a == 6);

	a = GetCrossAngle(6.1, 0.1, false);
	CHECK_APPROX(a, PI2 - 6.0);
}

TEST_CASE("Cross Time")
{
	const v2 center = { 0, 0 }; 
	const v2 pos    = { 3, 4 };
	const double r = 5;
	double ir = 0.1;
	double w = 2;

	double t = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, true); // @todo, pos/center omdraaien geeft GetPhi assert error!
	CHECK_APPROX_HIGH(t, 0.907 / 2);

	double t2 = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, false); // @todo, pos/center omdraaien geeft GetPhi assert error!
	CHECK_APPROX(t+t2, (PI2-(ir*2.0 / r))/ w);


	ir = 0.059;
	w = 100;
	
	t = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, true);
	t2 = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, false);
	CHECK_APPROX(t + t2, (PI2 - (ir * 2.0 / r)) / w);
}

#endif // TEST_H