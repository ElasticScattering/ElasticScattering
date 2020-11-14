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


double lifetime(const int quadrant, const int step, const double2 pos, BUFFER_ARGS);
double TraceOrbit(Particle* p, const Orbit* orbit, BUFFER_ARGS);
Orbit MakeOrbit(const double2 pos, const double phi, ScatteringParameters* sp);

inline Orbit MakeOrbit(const double2 pos, const double phi, ScatteringParameters* sp)
{
    const bool clockwise = sp->is_clockwise == 1;
    const bool incoherent = sp->is_incoherent == 1;
    const bool diag_regions = sp->is_diag_regions == 1;

    const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
    const double bound_angle = GetBoundAngle(phi, sp->alpha, clockwise);
    const double bound_phi = GetCrossAngle(phi, bound_angle, clockwise);

    const v2 vel = (double2)(cos(phi), sin(phi)) * sp->particle_speed;
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    const double2 center = GetCyclotronOrbitCenter(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    Orbit orbit(center, orbit_radius, clockwise, bound_time, bound_phi);

    return orbit;
}

inline double TraceOrbit(Particle* p, const Orbit* orbit, BUFFER_ARGS)
{
    double position_angle = GetPositionAngle(p->phi, orbit->clockwise); //@Todo, wat is dit...

    Intersection next_intersection;
    next_intersection.position = p->starting_position;
    next_intersection.entering_cell = p->cell_index;
    
    double lifetime = INF;
    while (1) {
        // Move to the next cell. @Todo, kopieert dit de referentie (conflict met volgende intersectie)?
        //Intersection entry_point = next_intersection;
        Intersection entry_point;
        entry_point.position = next_intersection.position;
        entry_point.entering_cell = next_intersection.entering_cell;
        entry_point.dphi = next_intersection.dphi;
        entry_point.incident_angle = next_intersection.incident_angle;

        // Find the next intersection point. If this orbit crosses the current
        // cell twice, we should limit valid intersections to those that happen before
        // leaving this cell for the next cell.
        double2 cell_pos = to_world(entry_point.entering_cell, sp->cells_per_row, sp->impurity_spawn_range);
        bool next_box_available = GetNextCell(orbit, cell_pos, entry_point, sp->cell_size, sp->cells_per_row, &next_intersection);
        double2 valid_phi_range = (double2)(entry_point.dphi, next_box_available ? next_intersection.dphi : position_angle); // ??

        // Try to find an intersection in the current cell.
        int impurity_start = (entry_point.entering_cell - 1 < 0) ? 0 : cell_indices[entry_point.entering_cell-1];
        int impurity_end = cell_indices[entry_point.entering_cell];

        for (int i = impurity_start; i < impurity_end; i++) {
            double2 impurity = impurities[i];
            if (CirclesCross(orbit, impurity, sp->impurity_radius)) {
                double t = GetFirstCrossTime(orbit, p->starting_position, impurity, sp->impurity_radius, sp->angular_speed, valid_phi_range);
                lifetime = (t < lifetime) ? t : lifetime;
            }
        }

        // Exit early if an intersection was found, or if no valid 
        // intersection can occur anymore, because the orbit doesn't
        // cross any more cells, or if the particle's lifetime has
        // ended.
        if (lifetime < INF || !next_box_available || orbit->bound_phi < GetCrossAngle(p->phi, next_intersection.dphi, orbit->clockwise))
            break;
    }

    return lifetime;
}

inline double lifetime(const int quadrant, const int step, const double2 pos, BUFFER_ARGS)
{
    const double phi = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;
    const Orbit orbit = MakeOrbit(pos, phi, sp);

    Particle p;
    p.starting_position = pos;
    p.phi = phi;
    p.cell_index = get_cell_index(pos, sp->impurity_spawn_range, sp->max_expected_impurities_in_cell);

    const double max_lifetime = min(sp->default_max_lifetime, orbit.bound_time);

    int impurity_start = (p.cell_index - 1 < 0) ? 0 : cell_indices[p.cell_index - 1];
    int impurity_end = cell_indices[p.cell_index];
    for (int i = impurity_start; i < impurity_end; i++)
    {
        if (InsideImpurity(pos, impurities[i], sp->impurity_radius))
            return 0;
    }

    double lt = TraceOrbit(&p, &orbit, sp, impurities, cell_indices);
    return lt;
}
