#include "src/sim/es/lifetime.h"
#include "src/sim/es/settings.h"
#include "src/sim/es/particle_metrics.h"

#include "src/sim/cl/cl_macros.h"

kernel void 
lifetime(constant SimulationSettings* ss, // Settings used for generic simulation settings, scaling etc.
		 constant ParticleSettings* ps,   // Settings used for creating the particle and orbit.
		 constant ImpuritySettings* is,   // Settings used for calculating intersections and grid movement.
		 constant double2* impurities, 
		 constant int* imp_index,
		 global double* lifetimes,
		 global ParticleMetrics* metrics)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
    
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;


	const int q = (int)(v / ss->particles_per_quadrant);
	const int p =       v % ss->particles_per_quadrant;

    double2 pos = (double2)(i, j) * ss->distance_between_positions + ss->small_offset;
	int idx = GET_INDEX(i, j, v);
	/*
#ifdef DEBUGGING
	if (idx < 100) {
		printf("Particle:\n");
		printf("i: %i j: %i v: %i\n", i, j, v);
		printf("x: %.8e y: %.8e phi: %.4e\n", pos.x, pos.y, GET_PHI(v));

		printf("\nParticle settings:\n");
		printf("is_clockwise: %i\n", ps->is_clockwise);
		printf("is_coherent: %i\n", ps->is_coherent);
		printf("particle_speed: %e\n", ps->particle_speed);
		printf("angular_speed: %e\n", ps->angular_speed);
		printf("alpha: %e\n", ps->alpha);
		printf("phi_start: %e\n", ps->phi_start);
		printf("phi_step_size: %e\n", ps->phi_step_size);

		printf("\nImpurity settings:\n");
		printf("spawn_region_start: %e\n", is->spawn_region_start);
		printf("spawn_region_size: %e\n", is->spawn_region_size);
		printf("cell_size: %e\n", is->cell_size);
		printf("cells_per_row: %i\n", is->cells_per_row);
		printf("impurity_radius: %e\n", is->impurity_radius);

		printf("\nSimulation settings:\n");
		printf("particles_per_quadrant: %i\n", ss->particles_per_quadrant);
		printf("particles_per_position: %i\n", ss->particles_per_position);
		printf("particles_per_row: %i\n", ss->particles_per_row);
		printf("total_particles: %i\n", ss->total_particles);

		printf("positions_per_row: %i\n", ss->positions_per_row);
		printf("total_positions: %i\n", ss->total_positions);

		printf("region_size: %e\n", ss->region_size);
		printf("region_extended_area: %e\n", ss->region_extended_area);
		printf("distance_between_positions: %e\n", ss->distance_between_positions);
		printf("small_offset: x: %e y: %e\n", ss->small_offset.x, ss->small_offset.y);

		printf("\nImpurities\n");

		for (int imp_idx = 0; imp_idx < 20; imp_idx++)
			printf("x: %.8e y: %.8e\n", impurities[imp_idx].x, impurities[imp_idx].y);
	}
#endif
*/

	Particle particle = CreateParticle(q, p, pos, ps);
	lifetimes[GET_INDEX(i, j, v)] = TraceOrbit(&particle, is, impurities, imp_index, metrics);
}
