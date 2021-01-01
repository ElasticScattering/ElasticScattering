// Oude variant, met alle phi waardes van één quadrant in één WI. Dit betekend dat er op ruime sprongen naar lifetimes geschreven wordt.
// Nadeel is dat phi/quadrant berekend moet worden, en ook de positie.

kernel void 
quadrant_lifetime(constant ParticleSettings* particle_settings,     // Settings used for creating the particle and orbit.
				  constant ImpuritySettings* impurity_settings,     // Settings used for calculating intersections and grid movement.
				  constant SimulationSettings* ss, // Settings used for generic simulation settings, scaling etc.
				  constant double2* impurities, 
				  constant int* imp_index,
				  global double* lifetimes
				  global Metrics* metrics)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int q = get_global_id(2);
    
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

    const double2 pos = ((double2)(i, j) * ss->distance_between_particles) + ss->small_offset;
	unsigned int base_idx = GET_INDEX(i, j, q);

	for (int p = 0; p < ss->integrand_steps; p++)
	{
		Particle particle = CreateParticle(q, p, pos, &particle_settings);
		lifetimes[base_idx + p] = TraceOrbit(&particle, &impurity_settings, impurities, imp_index, metrics);
	}
}

kernel void 
quadrant_apply_sigma_component(constant double* lifetimes, constant SimulationSettings* ss, constant ParticleSettings* ps, int mode, global double* sigma_component)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int q = get_global_id(2);
	
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

	//double tau = (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);
	
	unsigned int base_idx = GET_INDEX(i, j, q);

	for (int p = 0; p < sp->integrand_steps; p++)
	{
		double f = (mode == MODE_SIGMA_XX) cos(p) : sin(p); 
		sigma_component[base_idx + p] = GetSigma(lifetimes[base_idx + p], GET_PHI(q, p), tau, ss->signed_angular_speed);
	}
}


kernel void 
apply_simpson_weights_particles(global double* particle_values, global double* particle_results)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || y >= limit) {
		return;
	}

	const int idx = GET_PARTICLE_INDEX(i, j);
	particle_results[idx] = particle_values[idx] * SimpsonWeight2D(x, y, row_size-1);
}



kernel void 
integrate_to_particle(global double* values, double values_per_quadrant, double integrand_factor, global double* particle_results)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int row_size = get_global_size(0);
    
	int limit = row_size-1;
	if (i >= limit || y >= limit) {
		return;
	}

	double total = 0;
	for (int q = 0; q < 4; q++)
	{
		unsigned int base_idx = GET_INDEX(i, j, q);
	
		for (int p = 0; p < values_per_quadrant; p++)
			total += values[base_idx + p] * SimpsonWeight(p, values_per_quadrant);
	}
	
	particle_results[GET_PARTICLE_INDEX(i, j)] = total * integrand_factor;
}

