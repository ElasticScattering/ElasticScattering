__kernel void scatter(__global double2 *imps, 
                               double   imp_radius,
                      __global bool    *alive) 
{
    const int idx = get_global_id(0);

    for (int i = 0; i < 1000; i++) {
        alive[idx] = true;
    }
}