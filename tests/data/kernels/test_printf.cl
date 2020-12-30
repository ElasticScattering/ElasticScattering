__kernel void test_printf() {
	int id = get_global_id(0);
	printf("GPU id: %d", id);
}
