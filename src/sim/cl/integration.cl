/* -------------------------------------------------------------------------
	This code is part of ElasticScattering.
	
	Copyright(C) 2022 Elastic Scattering developers

	This program is free software : you can redistribute it and /or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#include "src/sim/es/util.h"
#include "src/sim/cl/cl_macros.h"
#include "src/sim/es/constants.h"
#include "src/sim/es/settings.h"

kernel void 
apply_max_lifetime(constant double* raw_lifetimes, const double default_max_lifetime, global double* lifetimes)
{
	int idx = get_global_id(0);

	double lt = raw_lifetimes[idx];
	lifetimes[idx] = min(lt, default_max_lifetime);
}

kernel void 
apply_sigma_component(constant double* lifetimes, 
					  constant SimulationSettings* ss, 
					  constant ParticleSettings* ps, 
					  const double tau, 
					  const int mode, 
					  global double* sigma_lifetimes)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
	
	if (i >= ss->positions_per_row || j >= ss->positions_per_row)
		return;

	double quadrant  = floor(v / (double)ss->particles_per_quadrant);
	double phi_index = (double)(v % ss->particles_per_quadrant); 
	double phi       = ps->phi_start + quadrant * HALF_PI + phi_index * ps->phi_step_size;

	int idx  = GET_INDEX(i, j, v);
	double f = (mode == MODE_SIGMA_XX) ? cos(phi) : sin(phi);
	sigma_lifetimes[idx] = f * GetSigma(lifetimes[idx], phi, tau, ss->signed_angular_speed);
}


kernel void 
apply_simpson_weights(global double* sigma_lifetimes, constant SimulationSettings* ss)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	int v = get_global_id(2);
    
	if (i >= ss->positions_per_row || j >= ss->positions_per_row)
		return;

	sigma_lifetimes[GET_INDEX(i, j, v)] *= (SimpsonWeight2D(i, j, ss->positions_per_row) * SimpsonWeight(v % ss->particles_per_quadrant, ss->particles_per_quadrant));
}

kernel void 
integrate_to_position(global double* values, constant SimulationSettings* ss, global double* position_values)
{
	int i = get_global_id(0);
    int j = get_global_id(1);
	
	if (i >= ss->positions_per_row || j >= ss->positions_per_row)
		return;

	int row_size = get_global_size(0);
	int base_idx = j * row_size * ss->particles_per_position + i * ss->particles_per_position;
	int pos_idx  = j * row_size + i;
	
	double total = 0;
	for (int q = 0; q < 4; q++)
		for (int p = 0; p < ss->particles_per_quadrant; p++)
			total += values[base_idx++] * SimpsonWeight(p, ss->particles_per_quadrant);
	
	position_values[pos_idx] = total * ss->phi_integrand_factor;
}

kernel void 
sum(global double* A, local double* local_sums, global double* B)
{
	uint id = get_global_id(0);
	uint local_id = get_local_id(0);

	local_sums[local_id] = A[id] + A[id + get_global_size(0)]; // Global size must be half the length of A.
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
