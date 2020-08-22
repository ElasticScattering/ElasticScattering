#include <CL/cl.h>

typedef struct
{
	cl_int dim;
	cl_int particle_count;
	cl_int impurity_count;

	cl_double region_size;
	cl_double particle_max_lifetime;
	cl_double particle_speed;      // v
	cl_double particle_mass;       // m
	cl_double impurity_radius;     // r
	cl_double impurity_radius_sq;  // r^2
	cl_double tau;

	cl_double alpha;
	cl_double phi;
	cl_double magnetic_field;      // B
	cl_double angular_speed;       // w
} SimulationParameters;