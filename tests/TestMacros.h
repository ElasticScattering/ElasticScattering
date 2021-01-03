#pragma once

#define EPSILON_RELATIVE 1e-5
#define EPSILON 1e-6
#define EPSILON_LOW_PRECISION 0.01
#define EPSILON_HIGH 1e-13

#define REQUIRE_ALMOST(a, b)  { REQUIRE(abs((a)-(b)) < EPSILON_HIGH); }

#define CHECK_ALMOST(a, b)  { CHECK(abs((a)-(b)) < EPSILON_HIGH); }
#define CHECK_ALMOST_MSG(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON_HIGH, p_msg); }

#define CHECK_APPROX(a, b)  { CHECK(abs((a)-(b)) < EPSILON); }
#define CHECK_APPROX_MSG(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON, p_msg); }
#define CHECK_APPROX_LOW(a, b)  { CHECK(abs((a)-(b)) < EPSILON_LOW_PRECISION); }

#define CHECK_RELATIVE_MF(a, b, mf) { CHECK(abs(((a)-(b)) / (b)) < (1.0 / mf) * EPSILON_RELATIVE ); }
#define CHECK_RELATIVE(a, b) { CHECK(abs(((a)-(b)) / (b)) < EPSILON_RELATIVE ); }