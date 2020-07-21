kernel void doubled(global int* A) {
    const int idx = get_global_id(0);
    A[idx] = A[idx] * 10;
}