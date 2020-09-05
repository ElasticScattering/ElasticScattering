#include "src/escl/constants.h"

__kernel void sum(__global double* A, __global double* B, __local double* local_sums)
{
	uint id = get_global_id(0);
	uint local_id = get_local_id(0);

	local_sums[local_id] = A[id] + A[id + get_global_size(0)]; // Global size should be half the length of A
	barrier(CLK_LOCAL_MEM_FENCE);
	for (int stride = get_local_size(0) / 2; stride > 1; stride /= 2) {
		if (local_id < stride) {
			local_sums[local_id] += local_sums[local_id + stride];
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (local_id == 0) {
		B[get_group_id(0)] = local_sums[0] + local_sums[1];
	}
}

/*
__kernel void sum2(__global double* data, __global double* output, __local double* partial_sums)
{
	int lid = get_local_id(0);
	int group_size = get_local_size(0);

	partial_sums[lid] = data[get_global_id(0)];
	barrier(CLK_LOCAL_MEM_FENCE);

	for (int i = group_size / 2; i > 0; i >>= 1) {
		if (lid < i) {
			partial_sums[lid] += partial_sums[lid + i];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (lid == 0) {
		output[get_group_id(0)] = partial_sums[0];
	}
}
*/

