#pragma once

#include "lifetime.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "impurity_grid.h"

#include "windows.h"


double SingleLifetime(const Particle* p, const double2 impurity, const double impurity_radius, const double angular_speed)
{
    double2 d = p->starting_position - impurity;
    if ((impurity_radius * impurity_radius) > dot(d, d))
        return 0;

    if (CirclesCross(p->orbit, impurity, impurity_radius))
    {
        return GetFirstCrossTime(p->starting_position, p->orbit, impurity, impurity_radius, angular_speed);
    }

    return DBL_MAX;
}

double TraceOrbit(Particle* p, BUFFER_ARGS)
{
    Intersection entry_point;
    entry_point.position = p->starting_position;
    entry_point.dphi = 0; // ????

    // Early exit als laatste intersectie een grotere hoek heeft dan bound_angle.

    double lifetime = DBL_MAX;
    bool hit = false;
    while (!hit) {
        int impurity_start = cell_indices[p->cell_index];
        int impurity_end = cell_indices[p->cell_index + 1]; // null?

        for (int i = impurity_start; i < impurity_end; i++) {
            double t = SingleLifetime(p, impurities[i], sp->impurity_radius, sp->angular_speed);

            lifetime = (t < lifetime) ? t : lifetime;
            hit = (lifetime < DBL_MAX);
        }

        if (!hit) {
            double2 cell_pos = to_world(p->cell_index, sp->cells_per_row, sp->impurity_spawn_range);

            int next_cell;
            Intersection next_intersection;
            bool success = GetNextCell(p->orbit, p->cell_index, cell_pos, entry_point, sp->cell_size, sp->cells_per_row, &next_cell, &next_intersection);

            if (!success) {
                // The End?
                break;
            }

            p->cell_index = next_cell;
            entry_point = next_intersection;
        }
    }

    return lifetime;
}

double lifetime(const int quadrant, const int step, const double2 pos, BUFFER_ARGS)
{
    const double phi = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;
    const Orbit orbit = MakeOrbit(pos, phi, sp);
   
    Particle p;
    p.starting_position = pos;
    p.orbit = orbit;
    p.phi = phi;
    p.cell_index = get_cell_index(pos, sp->impurity_spawn_range, sp->max_expected_impurities_in_cell);

    const double max_lifetime = min(sp->default_max_lifetime, orbit.bound_time);

    // @Todo, Check if we start inside an impurity.
    int impurity_start = cell_indices[p.cell_index];
    // ....

    double lt = TraceOrbit(&p, sp, impurities, cell_indices);
    return lt;
}


Orbit MakeOrbit(const double2 pos, const double phi, ScatteringParameters* sp)
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
