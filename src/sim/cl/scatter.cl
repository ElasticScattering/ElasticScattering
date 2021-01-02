#include "src/sim/es/lifetime.h"
#include "src/sim/es/settings.h"
#include "src/Metrics.h"

#include "src/sim/cl/cl_macros.h"

kernel void 
lifetime(constant SimulationSettings* ss, // Settings used for generic simulation settings, scaling etc.
		 constant ParticleSettings* ps,   // Settings used for creating the particle and orbit.
		 constant ImpuritySettings* is,   // Settings used for calculating intersections and grid movement.
		 constant double2* impurities, 
		 constant int* imp_index,
		 global double* lifetimes,
		 global Metrics* metrics)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
    
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

	const int q = (int)(v / ss->values_per_quadrant);
	const int p =       v % ss->values_per_quadrant;
    const double2 pos = ((double2)(i, j) * ss->distance_between_particles) + ss->small_offset;
	
	Particle particle = CreateParticle(q, p, pos, ps);
	lifetimes[GET_INDEX(i, j, v)] = TraceOrbit(&particle, is, impurities, imp_index, metrics);
}
