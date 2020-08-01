__kernel void scatter(__global double2 *imps, 
                               double   imp_radius) 
{
    const int idx = get_global_id(0);

    for (int i = 0; i < imp_radius; i++) {
    }
}