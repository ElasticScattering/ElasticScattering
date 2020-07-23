kernel void doubled(global int* A) {
    const int idx = get_global_id(0);
    A[idx] = A[idx] * 10;
}


kernel void double_precision(global double* A) {
    const int idx = get_global_id(0) % 1000;
    A[idx] = A[idx] * 10.03;
}