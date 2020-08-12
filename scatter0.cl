__constant double PI = 3.141592653589793238463;
__constant double PI2 = 6.283185307179586;

#define GLINTEROP

#define ROW_SIZE 1000


typedef struct Simulation {
    double region_size;
	int particle_count;
	int particle_row_count;
    double impurity_radius;     // r
	double impurity_radius_sq;  // r^2
    double magnetic_field;      // B
    double _max_lifetime;
	double speed;      // v
	double mass;       // m
	double tau;
    double alpha;
} Simulation;

double lifetime0(double tau, double2 pos, double phi, double speed, int impurity_count, double imp_radius, __global double2 *imps)
{
    double impurity_radius_sq = imp_radius * imp_radius;
    double lifetime = tau;
      
    double2 unit = { cos(phi), sin(phi) };
    double2 vel = unit * speed;
    
    for (int i = 0; i < impurity_count; i++) {
        double2 imp_pos = imps[i];
        double inner = (imp_pos.x - pos.x) * unit.x + (imp_pos.y - pos.y) * unit.y;
        double2 projected = pos + unit * inner;
        
        double2 d = projected - imp_pos;
        double diff = impurity_radius_sq - dot(d, d);
        if (diff < 0.0) {
            continue;
        }

        double L = sqrt(diff);

        double2 time_taken;
        if (vel.x != 0) {
            time_taken.x = -((projected.x - L * unit.x) - pos.x) / vel.x;
            time_taken.y = -((projected.x + L * unit.x) - pos.x) / vel.x;
        }
        else {
            time_taken.x = -((projected.y - L * unit.y) - pos.y) / vel.y;
            time_taken.y = -((projected.y + L * unit.y) - pos.y) / vel.y;
        }
        

        if ((time_taken.x * time_taken.y) < 0) {
            lifetime = 0;
            break;        
		}

        if (time_taken.x > 0 && time_taken.x < lifetime) {
            lifetime = time_taken.x;
        }
        if (time_taken.y > 0 && time_taken.y < lifetime) {
            lifetime = time_taken.y;
        }
    }

    return lifetime;
}
__kernel void lifetime(double region_size, double speed, double tau, double phi, int impurity_count, double imp_radius, __global double2 *imps,
#ifdef GLINTEROP
                       __write_only image2d_t screen)
#else
                       __global double *lifetimes) 
#endif
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    int row_size = get_global_size(0)-1;
    double2 pos = (double2)(region_size * x, region_size * y) / row_size;

    double lifetime = lifetime0(tau, pos, phi, speed, impurity_count, imp_radius, imps);

#ifdef GLINTEROP
    float k = (float)(lifetime / tau);
    write_imagef(screen, (int2)(x, y), (float4)(k,k,k,1.0f));
#else
    lifetimes[y * (row_size+1) + x] = lifetime;
#endif
}