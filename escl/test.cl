__kernel void many_sqrt(__global double *a, int n, __global double *b) {
	int id = get_global_id(0);
	for (int i = 0; i < n; i++) {
		double w = fabs(sin((double)(i)));
	    b[id] += a[id] * sqrt(w);
	}
}