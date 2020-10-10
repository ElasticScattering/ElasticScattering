#pragma once

#include "doctest.h"
#include "TestMacros.h"
#include "src/ElasticScattering.h"
#include "src/ParametersFactory.h"


/* TODO: Get full results not 1.231e07
TEST_CASE("Compare Sigma XX/XY results to verified results") {
    ScatteringParameters sp = ParametersFactory::GenerateNoImpurities();

    auto es = new GPUElasticScattering();




    for (int i = 0; i < 50; i++) {
        sp.magnetic_field = i;
        sp.mode = MODE_SIGMA_XX;
        double result = es->Compute(sp);

        CHECK(result == result);

        sp.mode = MODE_SIGMA_XY;

        double result2 = es->Compute(sp);

        CHECK(result2 == result2);
    }
}
*/

TEST_CASE("Compare Sigma XX to formula (no impurities)") {
	ScatteringParameters sp = ParametersFactory::GenerateNoImpurities();
    sp.mode = MODE_SIGMA_XX;

	auto e = new CPUElasticScattering;
    double result = e->Compute(sp);

	double kf = M * sp.particle_speed / HBAR;
	double n = (kf * kf) / (PI2 * C1);
	double formula = E * E * n * sp.tau / M;

	CHECK_RELATIVE(result, formula);
}
