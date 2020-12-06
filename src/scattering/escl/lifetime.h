#pragma once

#include "ScatteringParameters.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "impurity_grid.h"

#ifndef DEVICE_PROGRAM
    #include <vector>
    #include "windows.h"
#endif

/// <summary>
/// Returns the time at which the orbit collides with the first impurity, or INF if no collision happened. 
/// </summary>
ESCL_INLINE double TraceOrbit(const Particle* const p, const Orbit* const orbit, BUFFER_ARGS)
{
    double position_angle = GetPositionAngle(p->phi, orbit->clockwise); //@Todo, wat is dit...

    double lifetime = INF;
    Intersection next_intersection;
    next_intersection.position = p->starting_position;
    next_intersection.entering_cell = p->starting_cell;
    next_intersection.dphi = PI2;
    next_intersection.incident_angle = p->phi;
    
    while (1) {
        // Move to the next cell.
        Intersection entry_point = next_intersection;

        // Find the next intersection point. 
        bool next_box_available = GetNextCell(orbit, entry_point, sp->cell_size, sp->cells_per_row, sp->impurity_spawn_range, &next_intersection);
        
        // Not all impurites in a cell should be tested, because an orbit can 
        // leave this cell and then come back in later. Only the impurities 
        // that the current orbit segment can intersect with should be tested.
        double2 valid_phi_range = (double2)(entry_point.dphi, next_box_available ? next_intersection.dphi : position_angle);

        // Get the impurities in this cell from the index.
        int cell_idx = get_index(entry_point.entering_cell, sp->cells_per_row);
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
        if (lifetime < INF || !next_box_available || (sp->is_incoherent && orbit->bound_phi < GetCrossAngle(p->phi, next_intersection.dphi, orbit->clockwise)))
            break;
    }
    
    // @Todo, Hou bij of particle ergens op gebotst heeft (default_max_lifetime).
    
    return min(lifetime, orbit->bound_time);
}

ESCL_INLINE double lifetime(const int quadrant, const int step, const double2 pos, BUFFER_ARGS)
{
    Particle p;
    p.starting_position = pos;
    p.phi = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;
    p.starting_cell = get_cell(pos.x, pos.y, sp->impurity_spawn_range, sp->cells_per_row);

    int particle_cell_index = get_index(p.starting_cell, sp->cells_per_row);
    int impurity_start = (p.starting_cell == 0) ? 0 : cell_indices[particle_cell_index - 1];
    int impurity_end = cell_indices[particle_cell_index];
    for (int i = impurity_start; i < impurity_end; i++)
    {
        if (InsideImpurity(pos, impurities[i], sp->impurity_radius))
            return 0;
    }

    const bool clockwise = sp->is_clockwise == 1;
    const bool incoherent = sp->is_incoherent == 1;

    const double bound_time = GetBoundTime(p.phi, sp->alpha, sp->angular_speed, incoherent, clockwise, false);
    const double bound_angle = GetBoundAngle(p.phi, sp->alpha, clockwise);
    const double bound_phi = sp->is_incoherent ? GetCrossAngle(p.phi, bound_angle, clockwise) : INF;

    const double2 dir = { cos(p.phi), sin(p.phi) };
    const double2 vel = dir * sp->particle_speed;
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    const double2 center = GetCyclotronOrbitCenter(p.starting_position, vel, orbit_radius, sp->particle_speed, clockwise);

    Orbit orbit(center, orbit_radius, clockwise, bound_time, bound_phi);

    double lt = TraceOrbit(&p, &orbit, sp, impurities, cell_indices);
    return lt;
}
