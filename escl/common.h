#ifndef CL_COMMON_H
#define CL_COMMON_H

__constant double PI = 3.141592653589793238463;
__constant double PI2 = 6.283185307179586;

struct Circle
{
	double2 pos;
	double radius;
};


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
		B[get_group_id(0)] = local_sums[0] +  local_sums[1];
	}
}


__kernel void to_texture(__global double* lifetimes, __write_only image2d_t screen)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	int row_size = get_global_size(0);

	float k = (float)(lifetimes[y * row_size + x] / 1e-12); // @Todo, wat is hier de max lifetime?
	write_imagef(screen, (int2)(x, y), (float4)(k, k, k, 1.0f));
}


__kernel void sigma_xx(__global double* lifetimes)
{

}
#endif // CL_COMMON_H