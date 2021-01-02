#include "./tests/sim/cl/kernels/test_cpp.h"

__kernel void test_cpp(__global int *B) {
	int id = get_global_id(0);
	int result = id;
	increment_var2(&result);
	B[id] = result;
}
