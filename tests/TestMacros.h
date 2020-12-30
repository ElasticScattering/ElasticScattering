#pragma once

#define EPSILON 0.000001
#define EPSILON_LOW_PRECISION 0.01
#define EPSILON_HIGH 1e-13

#define REQUIRE_ALMOST(a, b)  { REQUIRE(abs((a)-(b)) < EPSILON_HIGH); }

#define CHECK_ALMOST(a, b)  { CHECK(abs((a)-(b)) < EPSILON_HIGH); }
#define CHECK_ALMOST_MSG(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON_HIGH, p_msg); }
#define CHECK_APPROX(a, b)  { CHECK(abs((a)-(b)) < EPSILON); }
#define CHECK_APPROX_MSG(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON, p_msg); }
#define CHECK_APPROX_LOW(a, b)  { CHECK(abs((a)-(b)) < EPSILON_LOW_PRECISION); }
#define CHECK_RELATIVE(a, b) { CHECK(abs(((a)-(b)) / (b)) < EPSILON); }

#define CHECK_CPU_GPU_ALMOST(p_msg)                      \
	e->Compute(sp, cpu_result);							 \
	e2->Compute(sp, gpu_result);						 \
	CHECK_ALMOST(cpu_result, gpu_result, p_msg);  

#define CHECK_CPU_GPU_APPROX(p_msg)                     \
	e->Compute(sp, cpu_result);                         \
	e2->Compute(sp, gpu_result);                        \
	CHECK_APPROX_MSG(cpu_result, gpu_result, p_msg);  
