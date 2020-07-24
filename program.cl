kernel void double_precision(global double2* A) {
    const int idx = get_global_id(0);
    A[idx].x = A[idx].x * 10.03;
}