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
double SingleLifetime(const Particle* p, const Orbit* orbit, const double2 impurity, const double impurity_radius, const double angular_speed, const double2 valid_range);
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

inline double SingleLifetime(const Particle* p, const Orbit* orbit, const double2 impurity, const double impurity_radius, const double angular_speed, const double2 valid_range)
{
    if (CirclesCross(orbit, impurity, impurity_radius))
    {
        return GetFirstCrossTime(orbit, p->starting_position, impurity, impurity_radius, angular_speed, valid_range);
    }

    return DBL_MAX;
}

inline double TraceOrbit(Particle* p, const Orbit* orbit, BUFFER_ARGS)
{
    Intersection entry_point;
    entry_point.position = p->starting_position;
    entry_point.dphi = 0; // ????

    // Early exit als laatste intersectie een grotere hoek heeft dan bound_angle.

    double lifetime = DBL_MAX;
    bool hit = false;
    while (!hit) {
        double2 cell_pos = to_world(p->cell_index, sp->cells_per_row, sp->impurity_spawn_range);
        
        int next_cell;
        Intersection next_cell_intersection;
        bool next_box_available = GetNextCell(orbit, p->cell_index, cell_pos, entry_point, sp->cell_size, sp->cells_per_row, &next_cell, &next_cell_intersection);

        double angle_max = next_box_available ? next_cell_intersection.dphi : GetPositionAngle(p->phi, orbit->clockwise);
        double2 valid_phi_range = (double2)(entry_point.dphi, angle_max);
        
        int impurity_start = (p->cell_index - 1 < 0) ? 0 : cell_indices[p->cell_index - 1];
        int impurity_end = cell_indices[p->cell_index];

        for (int i = impurity_start; i < impurity_end; i++) {
            double t = SingleLifetime(p, orbit, impurities[i], sp->impurity_radius, sp->angular_speed, valid_phi_range);

            lifetime = (t < lifetime) ? t : lifetime;
            hit = (lifetime < DBL_MAX);
        }

        if (!hit) {
            if (!next_box_available) {
                // The End?
                break;
            }

            // Move to the next cell.
            p->cell_index = next_cell;
            entry_point = next_cell_intersection;
        }
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
