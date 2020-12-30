#include "src/scattering/escl/lifetime.h"

//#define GET_INDEX(m_i, m_j, m_q, m_p) ((m_j * ss->values_per_row) + (m_i * ss->values_per_particle) + (m_q * ss->values_per_quadrant) + m_p)
#define GET_INDEX(m_i, m_j, m_v) ((m_j * ss->values_per_row) + (m_i * ss->values_per_particle) + m_v)
#define GET_PARTICLE_INDEX(i, j) (j * limit + i)

#define GET_PHI(m_q, m_p) (ps->phi_start + m_q * HALF_PI + m_p * ps->phi_step_size)



kernel void 
lifetime(constant SimulationSettings* ss, // Settings used for generic simulation settings, scaling etc.
		 constant ParticleSettings* ps,   // Settings used for creating the particle and orbit.
		 constant ImpuritySettings* is,   // Settings used for calculating intersections and grid movement.
		 constant double2* impurities, 
		 constant int* imp_index,
		 global double* lifetimes
		 global Metrics* metrics)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
    
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

	const int q = (int)(v / ss->integrand_steps);
	const int p =       v % ss->integrand_steps;
    const double2 pos = ((double2)(i, j) * ss->distance_between_particles) + ss->small_offset;
	
	Particle particle = CreateParticle(q, p, pos, &ps);
	lifetimes[GET_INDEX(i, j, v)] = TraceOrbit(&particle, &is, impurities, imp_index, metrics);
}

/* 

Oude variant, met alle phi waardes van één quadrant in één WI. Dit betekend dat er op ruime sprongen naar lifetimes geschreven wordt.
Nadeel is dat phi/quadrant berekend moet worden, en ook de positie.

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
*/

kernel void 
apply_sigma_component(constant double* lifetimes, constant SimulationSettings* ss, constant ParticleSettings* ps, int mode, global double* sigma_component)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
	
	if (i >= ss->particles_per_row || j >= ss->particles_per_row)
		return;

	unsigned int idx = GET_INDEX(i, j, v);
	const int q = (int)(v / ss->integrand_steps);
	const int p =       v % ss->integrand_steps;

	//double tau = (ps.is_coherent) ? ss.coherent_tau : HBAR / (KB * temperature);

	double f = (mode == MODE_SIGMA_XX) cos(p) : sin(p); 
	sigma_component[idx] = GetSigma(lifetimes[idx], GET_PHI(q, p), tau, ss->signed_angular_speed);
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
sum(global read_only double* A, global write_only double* B, local double* local_sums)
{
	uint id = get_global_id(0);
	uint local_id = get_local_id(0);

	local_sums[local_id] = A[id] + A[id + get_global_size(0)]; // Global size should be half the length of A
	barrier(CLK_LOCAL_MEM_FENCE);
	for (int stride = get_local_size(0) / 2; stride > 1; stride /= 2) {
		if (local_id < stride) {
			local_sums[local_id] += local_sums[local_id + stride];
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (local_id == 0) {
		B[get_group_id(0)] = local_sums[0] + local_sums[1];
	}
}
