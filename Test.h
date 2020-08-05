#ifndef TEST_H
#define TEST_H

//#include "Details.h"
#include <assert.h>

const double test_var = 1234.21;

void CyclotronTest() {
	v2 pos = { 1e-6, 1e-6 };

	v2 vel = { 7e5, 0 };
	double wc = E * 2.5 / M0;
	double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
	double radius = vf / wc;
	//auto c = GetCyclotronOrbit(pos, vel, radius, vf, true);
}




#endif // TEST_H