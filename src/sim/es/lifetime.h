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

#pragma once

#include "settings.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "cell_grid.h"

#include "particle_metrics.h"

#ifndef DEVICE_PROGRAM
    #include <vector>
    #include "windows.h"
#endif

ESCL_INLINE Particle CreateParticle(const int quadrant, const int phi_step_index, const double2 position, PARTICLE_SETTINGS)
{
    Particle p;
    p.position      = position;
    p.phi           = ps->phi_start + quadrant * HALF_PI + phi_step_index * ps->phi_step_size;
    p.angular_speed = ps->angular_speed;

    const double2 vel        = MAKE_DOUBLE2(cos(p.phi), sin(p.phi)) * ps->particle_speed;
    const double bound_angle = GetBoundAngle(p.phi, ps->alpha, ps->is_clockwise);

    Orbit orbit;
    orbit.clockwise      = ps->is_clockwise;
    orbit.radius         = ps->particle_speed / ps->angular_speed;
    orbit.radius_squared = orbit.radius * orbit.radius;
    orbit.center         = GetCyclotronOrbitCenter(p.position, vel, orbit.radius, ps->particle_speed, ps->is_clockwise);
    orbit.bound_time     = GetBoundTime(p.phi, ps->alpha, ps->angular_speed, ps->is_coherent, ps->is_clockwise, false);
    orbit.bound_phi      = ps->is_coherent ? INF : GetCrossAngle(p.phi, bound_angle, ps->is_clockwise);
    orbit.particle_angle = p.phi;

    p.orbit = orbit;

    return p;
}

/*
Returns the time at which the orbit collides with the first impurity, or INF 
if no collision happened.

Instead of looping through all impurities, we use the grid cells to only check
the impurities within cells that the particle will cross. Because we 
follow the particle's orbit, the first intersection we find this way, is also 
guaranteed to be the first actual intersection. This eliminates almost all
checks needed in simulations with a dense field of impurities.

There is an edge case where a particle might travel through a cell multiple times,
for example: once in the beginning and once at the end. If we find an intersection 
that actually happened in the end, we might conclude that that is the first 
intersection. However, an intersection in a subsequent cell could actually be earlier.
All this means is that if we do find an intersection, we must also verify that this 
was actually a valid intersection at that current stage. This can be checked by only 
consdering the impurities that could occur before the particle leaves the cell.
*/
ESCL_INLINE double TraceOrbit(const Particle* const p, IMPURITY_SETTINGS, BUFFER_ARGS)
{
    double lifetime = INF;

    Intersection next_intersection;
    next_intersection.position       = p->position;
    next_intersection.entering_cell  = get_cell(p->position, is->spawn_region_start, is->spawn_region_size, is->cells_per_row);
    next_intersection.dphi           = PI2;
    next_intersection.incident_angle = p->phi;
    
    int particle_cell_index = get_index(next_intersection.entering_cell, is->cells_per_row);

    int impurity_start = (particle_cell_index == 0) ? 0 : cell_indices[particle_cell_index - 1];	
    int impurity_end = cell_indices[particle_cell_index];

    for (int i = impurity_start; i < impurity_end; i++)
    {
        if (InsideImpurity(p->position, impurities[i], is->impurity_radius)) {
            METRIC_INC(particles_inside_impurity);
            return 0;
        }
    }

    while (1) {
        // Move to the next cell.
        Intersection entry_point = next_intersection;

        bool next_cell_available = GetNextCell(&p->orbit, p->phi, is, &entry_point, &next_intersection);
        double2 valid_phi_range = MAKE_DOUBLE2(entry_point.incident_angle, next_cell_available ? next_intersection.incident_angle : p->phi); // move to Intersection?

        // Use the grid index to get the impurities in the current cell.
        int cell_idx = get_index(entry_point.entering_cell, is->cells_per_row);
        int impurity_start = (cell_idx > 0) ? cell_indices[cell_idx - 1] : 0;
        int impurity_end = cell_indices[cell_idx];

        METRIC_INC(cells_passed);
        METRIC_ADD(impurities_tested, (impurity_end - impurity_start));

        // Test each impurity.
        for (int i = impurity_start; i < impurity_end; i++) {
            double2 impurity = impurities[i];
            if (CirclesCross(&p->orbit, impurity, is->impurity_radius)) {
                double t = GetFirstCrossTime(&p->orbit, impurity, is->impurity_radius, p->angular_speed, valid_phi_range);
                lifetime = (t < lifetime) ? t : lifetime;
            }
        }

        // Exit early if an intersection was found, or if no valid 
        // intersection can occur anymore or because the orbit doesn't
        // cross any more cells. 
        
        if (lifetime < INF) break;

        if (!next_cell_available)
        {
            METRIC_INC(particles_escaped);
            lifetime = 0;
            break;
        }

        // It would also be possible to stop coherent particles early if the particle's lifetime has ended.
        //bool bound_reached = !ps->is_coherent && orbit->bound_phi < GetCrossAngle(p->phi, next_intersection.dphi, orbit->clockwise);
    }

    if (lifetime > p->orbit.bound_time)
    {
        METRIC_INC(particles_at_bound);
        return p->orbit.bound_time;
    }
    return lifetime;
}
