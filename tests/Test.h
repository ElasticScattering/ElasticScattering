#ifndef TEST_H
#define TEST_H

#include "src/sim/Simulation.h"
#include "src/sim/es/lifetime.h"
#include "src/sim/es/util.h"
#include "src/utils/OpenCLUtils.h"
#include "src/utils/ParametersFactory.h"

#include <assert.h>
#include <random>
#include <limits>

#include "doctest.h"
#include "TestMacros.h"


/*
TEST_CASE("Different scatter modes")
{
	ScatteringParameters sp = ParametersFactory::GenerateDefault();
	
	auto e  = new SimulationCPU();
	auto e2 = new GPUSimulation();

	double cpu_result, gpu_result, diff;
	
	sp.mode = MODE_DIR_LIFETIME;
	CHECK_CPU_GPU_ALMOST("Direct Lifetime")

	sp.mode = MODE_PHI_LIFETIME;
	CHECK_CPU_GPU_ALMOST("Phi Lifetime")

	sp.mode = MODE_SIGMA_XX;
	CHECK_CPU_GPU_APPROX("Sigma XX Lifetime")

	sp.mode = MODE_SIGMA_XY;
	CHECK_CPU_GPU_APPROX("Sigma XY Lifetime")
}
*/




/*
TEST_CASE("Comparing kernel results on CPU and GPU")
{
	auto e  = new SimulationCPU();
	auto e2 = new GPUSimulation();

	ScatteringParameters sp;
	sp.region_size     = 1e-6;
	sp.dim             = 64;
	sp.particle_speed  = 7e5;
	sp.impurity_count  = 100;
	sp.impurity_radius = 1.5e-8;
	sp.alpha           = PI / 4.0;
	sp.phi             = 0;
	sp.magnetic_field  = 0;
	sp.tau             = 1e-12;
	sp.integrand_steps = 9;
	sp.is_clockwise    = 1;
	sp.region_extends  = sp.particle_speed * sp.tau;
	sp.is_incoherent = true;

	sp.mode            = MODE_DIR_LIFETIME;
	sp.impurity_seed   = 0;

	double cpu_result, gpu_result, diff;

	SUBCASE("Directional lifetime_old") {
		CHECK_CPU_GPU_ALMOST("Default parameters")

		sp.phi = -sp.alpha - 1e-10;
		CHECK_CPU_GPU_ALMOST("Different phi")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_ALMOST("More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_ALMOST("Larger impurities")

		sp.impurity_seed = 1;
		CHECK_CPU_GPU_ALMOST("Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_ALMOST("Magnetic field on")

		sp.is_clockwise = 0;
		CHECK_CPU_GPU_ALMOST("Clockwise off")
	}

	SUBCASE("Phi integrated lifetime_old") {
		sp.mode = MODE_PHI_LIFETIME;

		CHECK_CPU_GPU_APPROX("PHI - Default parameters")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_APPROX("PHI - More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_APPROX("PHI - Larger impurities")

		sp.impurity_seed = 2;
		CHECK_CPU_GPU_ALMOST("PHI - Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_APPROX("PHI - Magnetic field on")

		sp.is_clockwise = 0;
		CHECK_CPU_GPU_APPROX("PHI - Clockwise off")
	}
	

	SUBCASE("Sigma XX") {
		sp.mode = MODE_SIGMA_XX;

		CHECK_CPU_GPU_APPROX("SXX - Default parameters")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_APPROX("SXX - More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_APPROX("SXX - Larger impurities")

		sp.impurity_seed = 3;
		CHECK_CPU_GPU_ALMOST("SXX - Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_APPROX("SXX - Magnetic field on")

		sp.is_clockwise = 0;
		CHECK_CPU_GPU_APPROX("SXX - Clockwise off")
	}
	
	SUBCASE("Sigma XY") {
		sp.mode = MODE_SIGMA_XY;

		CHECK_CPU_GPU_APPROX("SXY - Default parameters")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_APPROX("SXY - More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_APPROX("SXY - Larger impurities")

		sp.impurity_seed = 3;
		CHECK_CPU_GPU_ALMOST("SXY - Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_APPROX("SXY - Magnetic field on")

		sp.is_clockwise = 0;
		CHECK_CPU_GPU_APPROX("SXY - Clockwise off")
	}
}
*/
#endif // TEST_H