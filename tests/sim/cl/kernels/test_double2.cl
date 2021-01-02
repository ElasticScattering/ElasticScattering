__kernel void test_double2(__global double *a) {
    int id = get_global_id(0);
    double2 x = (double2)(id, id) * 1.5;

	a[id] = x.x + x.y;
}