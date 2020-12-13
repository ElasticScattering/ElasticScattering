#pragma once

#include "ScatteringParameters.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "cell_grid.h"

#ifndef DEVICE_PROGRAM
    #include <vector>
    #include "windows.h"
#endif

/// <summary>
/// Returns the time at which the orbit collides with the first impurity, or INF if no collision happened.
/// </summary>
ESCL_INLINE double TraceOrbit(const Particle* const p, const Orbit* const orbit, BUFFER_ARGS)
{
    double lifetime = INF;

    Intersection next_intersection;
    next_intersection.position       = p->starting_position;
    next_intersection.entering_cell  = p->starting_cell;
    next_intersection.dphi           = PI2;
    next_intersection.incident_angle = p->phi;
    
    while (1) {
        // Move to the next cell.
        Intersection entry_point = next_intersection;

        /*
        Instead of looping through all impurities, we subdivide the grid 
        in a number of cells and only check impurities inside a single
        cell. 
        
        To ensure that the correct impurity is always the first intersection
        we encounter, we need to move through the cells following the 
        particle's orbit.
        
        Additionally, a particle might travel through a cell multiple times, which means 
        an impurity intersection in a later cell could happen before one in this 
        cell. To prevent this, only the impurities that happen before the particle
        will leave the cell are considered. 
        */
        bool next_cell_available = GetNextCell(orbit, p->phi, sp->cell_size, sp->cells_per_row, sp->impurity_spawn_range, &entry_point, &next_intersection);
        double2 valid_phi_range = v2(entry_point.dphi, next_cell_available ? next_intersection.dphi : p->position_angle);

        // Use the grid index to get the impurities in the current cell.
        int cell_idx       = get_index(entry_point.entering_cell, sp->cells_per_row);
        int impurity_start = (cell_idx > 0) ? cell_indices[cell_idx -1] : 0;
        int impurity_end   = cell_indices[cell_idx];

        // Test each impurity.
        for (int i = impurity_start; i < impurity_end; i++) {
            double2 impurity = impurities[i];
            if (CirclesCross(orbit, impurity, sp->impurity_radius)) {
                // @Todo, hoort starting_position hier? Voelt vreemd.
                double t = GetFirstCrossTime(orbit, p->starting_position, impurity, sp->impurity_radius, sp->angular_speed, valid_phi_range);
                lifetime = (t < lifetime) ? t : lifetime;
            }
        }

        // Exit early if an intersection was found, or if no valid 
        // intersection can occur anymore, because the orbit doesn't
        // cross any more cells, or if the particle's lifetime has
        // ended.
        if (lifetime < INF || !next_cell_available || (sp->is_incoherent && orbit->bound_phi < GetCrossAngle(p->phi, next_intersection.dphi, orbit->clockwise)))
            break;
    }
    
    return min(lifetime, orbit->bound_time);
}

ESCL_INLINE double lifetime(const int quadrant, const int step, const double2 pos, BUFFER_ARGS)
{
    const bool clockwise  = sp->is_clockwise == 1;
    const bool incoherent = sp->is_incoherent == 1;

    Particle p;
    p.starting_position = pos;
    p.phi               = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;
    p.starting_cell     = get_cell(pos, sp->impurity_spawn_range, sp->cells_per_row);
    p.position_angle    = GetPositionAngle(p.phi, clockwise);

    int particle_cell_index = get_index(p.starting_cell, sp->cells_per_row);
    int impurity_start      = (p.starting_cell == 0) ? 0 : cell_indices[particle_cell_index - 1];
    int impurity_end        = cell_indices[particle_cell_index];
    for (int i = impurity_start; i < impurity_end; i++)
    {
        if (InsideImpurity(pos, impurities[i], sp->impurity_radius))
            return 0;
    }

    const double bound_time = GetBoundTime(p.phi, sp->alpha, sp->angular_speed, incoherent, clockwise, false);
    const double bound_angle = GetBoundAngle(p.phi, sp->alpha, clockwise);
    const double bound_phi = sp->is_incoherent ? GetCrossAngle(p.phi, bound_angle, clockwise) : INF;

    const double2 dir = v2(cos(p.phi), sin(p.phi));
    const double2 vel = dir * sp->particle_speed;
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    const double2 center = GetCyclotronOrbitCenter(p.starting_position, vel, orbit_radius, sp->particle_speed, clockwise);

    Orbit orbit(center, orbit_radius, clockwise, bound_time, bound_phi);

    double lt = TraceOrbit(&p, &orbit, sp, impurities, cell_indices);
    return lt;
}
