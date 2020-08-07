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

typedef struct ParticleParameters {
	double phi;
	double angular_speed;       // w
} ParticleParameters;

__kernel void scatter0(double region_size,
                       double max_lifetime,
                       double speed,
                       double mass,
                       double imp_radius,
                       double tau,
                       double alpha,
                       double phi,
                       double magnetic_field,
                       double angular_speed,
                       int impurity_count,
                       __global double2 *imps,
#ifdef GLINTEROP
                       write_only image2d_t screen)
#else
                       __global double *lifetimes) 
#endif
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    double2 pos = {region_size * x / ROW_SIZE, region_size * y / ROW_SIZE};

    double2 unit = { cos(phi), sin(phi) };
    double2 vel = { speed * unit.x, speed * unit.y };

    double impurity_radius_sq = imp_radius * imp_radius;


    double lifetime = max_lifetime;

    for (int i = 0; i < impurity_count; i++) {
        double2 imp_pos = imps[i];
        double inner = (imp_pos.x - pos.x) * unit.x + (imp_pos.y - pos.y) * unit.y;
        double2 projected = { pos.x + inner * unit.x, pos.y + inner * unit.y };
        
        double d1 = projected.x - imp_pos.x;
        double d2 = projected.y - imp_pos.y;
        double a = d1*d1 + d2*d2;
        if (a > impurity_radius_sq) {
            continue;
        }

        double L = sqrt(impurity_radius_sq - a);

        double2 time_taken;
        if (vel.x != 0) {
            time_taken.x = -((projected.x - L * unit.x) - pos.x) / vel.x;
            time_taken.y = -((projected.x + L * unit.x) - pos.x) / vel.x;
        }
        else {
            time_taken.x = -((projected.y - L * unit.y) - pos.y) / vel.y;
            time_taken.y = -((projected.y + L * unit.y) - pos.y) / vel.y;
        }

        if ((time_taken.x * time_taken.y) < 0)
        {
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

#ifdef GLINTEROP
    double k = lifetime / max_lifetime;
    write_imagef(screen, (int2)(x, y), float4(k,k,k,1));
#else
    lifetimes[y * ROW_SIZE + x] = lifetime;
#endif
}