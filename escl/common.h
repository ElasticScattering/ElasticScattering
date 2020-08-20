#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "common_structs.h"

__constant double PI = 3.141592653589793238463;
__constant double PI2 = 6.283185307179586;
__constant double M0 = 9.109e-31;
__constant double E = 1.602e-19;
__constant double HBAR = 1.055e-34;
__constant double C = 1.15e-9;

__kernel void sum(__global double *A, __global double *B, __local double *local_sums)
{
	uint id       = get_global_id(0);
	uint local_id = get_local_id(0);
	
	local_sums[local_id] = A[id] + A[id+get_global_size(0)]; // Global size should be half the length of A
	barrier(CLK_LOCAL_MEM_FENCE);
	for (int stride = get_local_size(0)/2; stride > 1; stride/=2 ) {
		if (local_id < stride) {
			local_sums[local_id] += local_sums[local_id+stride];
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (local_id == 0) {
		B[get_group_id(0)] = local_sums[0] + local_sums[1];
	}
}

__kernel void to_texture(__global double* lifetimes, double tau, __write_only image2d_t screen)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	int row_size = get_global_size(0);

	float k = (float)(lifetimes[y * row_size + x]/tau);
	float4 c = (float4)(k, k, k, 1.0f);
	
	/*
	if (k > 1e-12)
		c = (float4)(0, 0, (k-1.5e-12)/5e-13, 1.0f);
	else
		c = (float4)(0, 0, 0, 1.0f);
	*/
	
	write_imagef(screen, (int2)(x, y), c);
}

/*
__kernel void struct_test(__global StructTest* st)
{
	st.test1 *= st.test2;
}
*/

#endif // CL_COMMON_H