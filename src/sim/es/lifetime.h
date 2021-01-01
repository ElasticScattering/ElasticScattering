#pragma once

#include "ScatteringParameters.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "cell_grid.h"

#include "src/Metrics.h"

#ifndef DEVICE_PROGRAM
    #include <vector>
    #include "windows.h"
#endif

ESCL_INLINE Particle CreateParticle(const int quadrant, const int step, const double2 pos, PARTICLE_SETTINGS)
{
    Particle p;
    p.starting_position = pos;
    p.phi               = ps->phi_start + quadrant * HALF_PI + step * ps->phi_step_size;
    p.angular_speed     = ps->angular_speed;

    const double2 vel        = MAKE_DOUBLE2(cos(p.phi), sin(p.phi)) * ps->particle_speed;
    const double bound_angle = GetBoundAngle(p.phi, ps->alpha, ps->is_clockwise);

    Orbit orbit;
    orbit.clockwise      = ps->is_clockwise;
    orbit.radius         = ps->particle_speed / ps->angular_speed;
    orbit.radius_squared = orbit.radius * orbit.radius;
    orbit.center         = GetCyclotronOrbitCenter(p.starting_position, vel, orbit.radius, ps->particle_speed, ps->is_clockwise);
    orbit.bound_time     = GetBoundTime(p.phi, ps->alpha, ps->angular_speed, ps->is_coherent, ps->is_clockwise, false);
    orbit.bound_phi      = ps->is_coherent ? INF : GetCrossAngle(p.phi, bound_angle, ps->is_clockwise);
    orbit.particle_angle = p.phi; //GetAngle(p.starting_position, &orbit); //@Refactor: name?

    p.orbit = orbit;

    return p;
}

/// <summary>
/// Returns the time at which the orbit collides with the first impurity, or INF if no collision happened.
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
/// </summary>
ESCL_INLINE double TraceOrbit(const Particle* const p, IMPURITY_SETTINGS, BUFFER_ARGS)
{
    double lifetime = INF;

    Intersection next_intersection;
    next_intersection.position       = p->starting_position;
    next_intersection.entering_cell  = get_cell(p->starting_position, settings->spawn_region_start, settings->spawn_region_size, settings->cells_per_row);
    next_intersection.dphi           = PI2;
    next_intersection.incident_angle = p->phi;
    
    int particle_cell_index = get_index(next_intersection.entering_cell, settings->cells_per_row);

    int impurity_start = (particle_cell_index == 0) ? 0 : cell_indices[particle_cell_index - 1];	
    int impurity_end = cell_indices[particle_cell_index];
    for (int i = impurity_start; i < impurity_end; i++)
    {
        if (InsideImpurity(p->starting_position, impurities[i], settings->impurity_radius)) {
            METRIC_INC(metrics->particles_inside_impurity);
            return 0;
        }
    }

    while (1) {
        // Move to the next cell.
        Intersection entry_point = next_intersection;

        bool next_cell_available = GetNextCell(&p->orbit, p->phi, settings, &entry_point, &next_intersection);
        double2 valid_phi_range = MAKE_DOUBLE2(entry_point.incident_angle, next_cell_available ? next_intersection.incident_angle : p->phi); // move to Intersection?

        // Use the grid index to get the impurities in the current cell.
        int cell_idx = get_index(entry_point.entering_cell, settings->cells_per_row);
        int impurity_start = (cell_idx > 0) ? cell_indices[cell_idx - 1] : 0;
        int impurity_end = cell_indices[cell_idx];

        METRIC_INC(metrics->cells_passed);
        METRIC_ADD(metrics->impurities_tested, (impurity_end - impurity_start));

        // Test each impurity.
        for (int i = impurity_start; i < impurity_end; i++) {
            double2 impurity = impurities[i];
            if (CirclesCross(&p->orbit, impurity, settings->impurity_radius)) {
                double t = GetFirstCrossTime(&p->orbit, impurity, settings->impurity_radius, p->angular_speed, valid_phi_range);
                lifetime = (t < lifetime) ? t : lifetime;
            }
        }

        // Exit early if an intersection was found, or if no valid 
        // intersection can occur anymore, because the orbit doesn't
        // cross any more cells, or if the particle's lifetime has
        // ended.
        
        //bool bound_reached = !ps->is_coherent && orbit->bound_phi < GetCrossAngle(p->phi, next_intersection.dphi, orbit->clockwise);
        if (lifetime < INF || !next_cell_available)
            break;
    }

    METRIC_INC_COND(lifetime > 1, metrics->particles_escaped);

    return min(lifetime, p->orbit.bound_time);
}
