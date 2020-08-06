#ifndef TEST_H
#define TEST_H

#include "Details.h"
#include <assert.h>

#define EPSILON 0.00001
#define EPSILONS 0.01
#define AssertEqual(a, b)  { assert(abs((a)-(b)) < EPSILON); }
#define AssertEqualS(a, b)  { assert(abs((a)-(b)) < EPSILONS); }

void CyclotronTest() {
	v2 pos = { 1e-6, 1e-6 };

	v2 vel = { 7e5, 0 };
	double wc = E * 2.5 / M0;
	double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
	double radius = vf / wc;
	auto c = GetCyclotronOrbit(pos, vel, radius, vf, true);
	AssertEqual(c.x, pos.x);
	AssertEqual(c.y, pos.y + radius);

	c = GetCyclotronOrbit(pos, vel, radius, vf, false);
	AssertEqual(c.x, pos.x);
	AssertEqual(c.y, pos.y - radius);
}

void BoundTimeDirectionTest()
{
	double phi = -0.25;
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(phi, w, alpha, false, true);
	AssertEqual(t, 0.25);

	t = GetBoundTime(phi, w, alpha, true, true);
	AssertEqual(t, 0.75);

	t = GetBoundTime(phi, w, alpha, true, false);
	AssertEqual(t, 0.25);

	t = GetBoundTime(phi, w, alpha, false, false);
	AssertEqual(t, 0.75);
}

void BoundTimeSectorsTest()
{
	double alpha = 1;
	double w = 0.5;
	
	double t = GetBoundTime(PI/2 - 0.25, w, alpha, false, true);
	AssertEqual(t, 0.25);

	t = GetBoundTime(PI - 0.25, w, alpha, false, true);
	AssertEqual(t, 0.25);

	t = GetBoundTime(3 * PI/2 - 0.25, w, alpha, false, true);
	AssertEqual(t, 0.25);
}

void BoundTimeCornerTest() 
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(0.5, w, alpha, false, true);
	AssertEqual(t, 1);

	t = GetBoundTime(0.5, w, alpha, false, false);
	AssertEqual(t, 0);
}

void CircleCrossTest()
{
	v2 p1 = { 1, 2 };
	v2 p2 = { 5, 5 };
	v2 p3 = { 1, 3 };
	assert(!CirclesCross(p1, 1, p2, 1));
	assert(!CirclesCross(p1, 1, p3, 5));
	assert(CirclesCross({ 0, 0 }, 1, { 0, 0.5 }, 0.6));
	assert(CirclesCross({ 0, 0 }, 1000, { 2, 0 }, 1001));
}

void CircleCrossPointsTest()
{
	// Symmetric
	auto ps = GetCrossPoints({ -1, 0 }, 1.5, { 1, 0 }, 1.5);
	v2 q1 = ps.first;
	v2 q2 = ps.second;
	AssertEqual(q1.x, 0);
	AssertEqual(abs(q1.y), sqrt(1.5*1.5 - 1));
	
	AssertEqual(q2.x, 0);
	AssertEqual(abs(q2.y), sqrt(1.5 * 1.5 - 1));

	// Somewhat symmetric
	ps = GetCrossPoints({ 0, 0 }, 1, { 1, 1 }, 1.5);
	q1 = ps.first;
	q2 = ps.second;
	AssertEqual(1, pow(q1.x, 2) + pow(q1.y, 2));
	
	AssertEqual(pow(1.5, 2), pow(q1.x -1, 2) + pow(q1.y-1, 2));
	
	AssertEqual(1, pow(q2.x, 2) + pow(q2.y, 2));
	AssertEqual(pow(1.5, 2), pow(q2.x - 1, 2) + pow(q2.y - 1, 2));

	// Asymmetric
	ps = GetCrossPoints({ 100, 0 }, 100, { -1, 0 }, 1.5);
	q1 = ps.first;
	q2 = ps.second;

	AssertEqual(pow(100, 2), pow(q1.x - 100, 2) + pow(q1.y, 2));
	AssertEqual(pow(1.5, 2), pow(q1.x +   1, 2) + pow(q1.y, 2));

	AssertEqual(pow(100, 2), pow(q2.x - 100, 2) + pow(q2.y, 2));
	AssertEqual(pow(1.5, 2), pow(q2.x +   1, 2) + pow(q2.y, 2));

	// Asymmetric, inside
	ps = GetCrossPoints({ 100, -0.5 }, 100, { 1, 0.5 }, 1.5);
	q1 = ps.first;
	q2 = ps.second;
	AssertEqual(pow(100, 2), pow(q1.x - 100, 2) + pow(q1.y + 0.5, 2));
	AssertEqual(pow(1.5, 2), pow(q1.x -   1, 2) + pow(q1.y - 0.5, 2));

	AssertEqual(pow(100, 2), pow(q2.x - 100, 2) + pow(q2.y + 0.5, 2));
	AssertEqual(pow(1.5, 2), pow(q2.x -   1, 2) + pow(q2.y - 0.5, 2));
}

void PhiTest()
{
	const v2 p1 = { 0, 0 };
	double phi = GetPhi({ 1, 0 }, p1, 1);
	AssertEqual(phi, 0);

	phi = GetPhi({ 0, 1 }, p1, 1);
	AssertEqual(phi, PI/2);

	phi = GetPhi({ -1, 0 }, p1, 1);
	AssertEqual(phi, PI);

	phi = GetPhi({ 0, -1 }, p1, 1);
	AssertEqual(phi, 3 * PI / 2);

	phi = GetPhi({ 0.99999, -0.0045 }, p1, 1);
	AssertEqualS(phi, PI2);
}

void CrossAngleTest()
{
	double a = GetCrossAngle(6.1, 0.1, true);
	AssertEqual(a, 6);
	
	a = GetCrossAngle(6.1, 0.1, false);
	AssertEqualS(a, 0.28);
}

void CrossTimeTest()
{
	const v2 center = { 3, 4 };
	const v2 pos    = { 0, 0 };
	double t = GetFirstCrossTime(pos, center, { 5, 0 }, 5, 0.1, 2, true); // @todo, pos/center omdraaien geeft GetPhi assert error!
	AssertEqualS(t, 0.907 / 2);

	//double t2 = GetFirstCrossTime(pos, center, { 5, 0 }, 5, 0.1, 2, false);
	//AssertEqualS(t+t2, PI);
}

void BoundTimeTest()
{
	BoundTimeDirectionTest();
	BoundTimeSectorsTest();
	BoundTimeCornerTest();
}

void RunAllTests()
{
	CyclotronTest();
	BoundTimeTest();
	CircleCrossTest();
	CircleCrossPointsTest();
	PhiTest();
	CrossAngleTest();
	CrossTimeTest();
}



#endif // TEST_H