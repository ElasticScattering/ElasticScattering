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

__kernel void add_integral_weights_2d(__global double* A)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	int row_size = get_global_size(0);

	bool is_padding = (x == (row_size - 1)) || (y == (row_size - 1));
	bool is_edge = (x == 0) || (x == (row_size - 2)) || (y == 0) || (y == (row_size - 2));

	double w = is_padding ? 0.0 : 1.0;
	if (!is_edge)
	{
		w = ((x % 2) == 0) ? 2.0 : 4.0;
		w *= ((y % 2) == 0) ? 2.0 : 4.0;
	}

	A[y * row_size + x] = w * A[y * row_size + x];
}

__kernel void to_texture(__global double* lifetimes, double tau, __write_only image2d_t screen)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	int row_size = get_global_size(0);
	float k = (float)(lifetimes[y * row_size + x]);

#if 1
	k /= tau;
	float4 c = (float4)(k, k, k, 1.0f);
#else
	float4 c;
	if (fabs(k - tau) < 0.00001 * tau)
		c = (float4)(0, 0, 1, 1.0f);
	else
		c = (float4)(k / tau, 0, 0, 1.0f);
#endif

	write_imagef(screen, (int2)(x, y), c);
}