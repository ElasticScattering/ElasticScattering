#pragma once

#include "ScatteringParameters.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "impurity_grid.h"
#include "orbit.h"

#include "windows.h"

inline double trace_grid(double2 pos, double phi, CellRange current_cell, BUFFER_ARGS)
{
    const v2 unit = { cos(phi), sin(phi) };
    const v2 vel = unit * sp->particle_speed;
    bool clockwise = (sp->is_clockwise == 1);

    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    const double2 center = GetCyclotronOrbitCenter(pos, vel, orbit_radius, sp->particle_speed, clockwise);
    const Orbit orbit(center, orbit_radius);

    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = 1112.0;
    bool hit = false;
    while (!hit) {
        current_cell = get_cell_range(pos, sp->impurity_spawn_range, sp->max_expected_impurities_in_cell);

        double low_left_x = sp->impurity_spawn_range.y - sp->impurity_spawn_range.x;

        double2 lower_left = { };

        int impurity_start = cell_indices[current_cell.start];
        int impurity_end = cell_indices[current_cell.end];

        for (int i = impurity_start; i < impurity_end; i++) {
            const double2 impurity = impurities[i];

            double2 d = pos - impurity;
            if (impurity_radius_sq > dot(d, d))
            {
                lifetime = 0;
                break;
            }

            if (CirclesCross(orbit, impurity, sp->impurity_radius))
            {
                double t = GetFirstCrossTime(pos, orbit, impurity, sp->impurity_radius, sp->angular_speed, clockwise);

                lifetime = (t < lifetime) ? t : lifetime;
                hit = true;
            }
        }

        if (!hit) {
            current_cell = GetNextCell();

        }
    }

    return lifetime;
}

inline double lifetime(const double max_lifetime, const double2 pos, const double phi, BUFFER_ARGS)
{
    CellRange current_cell = get_cell_range(pos, sp->impurity_spawn_range, sp->max_expected_impurities_in_cell);
    
    // @Todo, Check if we start inside an impurity.
    int impurity_start = cell_indices[current_cell.start];
    // ....

    // @Todo, bound intersect?
    
    //max_lifetime = 15.0 * sp->tau;  Move to SP

    double lt = trace_grid(pos, phi, current_cell, sp, impurities, cell_indices);

    return lt;
}

inline double sim_phi_lifetime(const double2 pos, int quadrant, int step, BUFFER_ARGS)
{
    bool clockwise = (sp->is_clockwise == 1);
    bool incoherent = (sp->is_incoherent == 1);
    bool diag_regions = (sp->is_diag_regions == 1);

    const double phi = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;

    const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
    return lifetime(min(15.0 * sp->tau, bound_time), pos, phi, sp, impurities, cell_indices);
}


//////// ! Deprecated
inline double lifetime_old(const double max_lifetime, const double2 pos, const double phi, const bool clockwise, BUFFER_ARGS)
{
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = { cos(phi), sin(phi) };
    vel = vel * sp->particle_speed;
    const double2 center = GetCyclotronOrbitCenter(pos, vel, orbit_radius, sp->particle_speed, clockwise);
    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime_old = max_lifetime;

    for (int i = 0; i < sp->impurity_count; i++) {
        const double2 impurity = impurities[i];

        double2 d = pos - impurity;
        if (impurity_radius_sq > dot(d, d))
        {
            lifetime_old = 0;
            break;
        }

        if (CirclesCross(center, orbit_radius, impurity, sp->impurity_radius))
        {
            double t = GetFirstCrossTime(center, pos, impurity, orbit_radius, sp->impurity_radius, sp->angular_speed, clockwise);

            lifetime_old = (t < lifetime_old) ? t : lifetime_old;
        }
    }

    return lifetime_old;
}