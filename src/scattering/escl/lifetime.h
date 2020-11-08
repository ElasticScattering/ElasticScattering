#pragma once

#include "ScatteringParameters.h"

#include "device_macros.h"
#include "details.h"
#include "constants.h"
#include "impurity_grid.h"

#include "windows.h"

inline double TraceOrbit(const double2 pos, double phi, const Orbit orbit, int current_cell, BUFFER_ARGS)
{
    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double2 entry_point = pos;
    
    double lifetime = DBL_MAX;
    bool hit = false;
    while (!hit) {
        int impurity_start = cell_indices[current_cell];
        int impurity_end   = cell_indices[current_cell+1];

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
                double t = GetFirstCrossTime(pos, orbit, impurity, sp->impurity_radius, sp->angular_speed);

                lifetime = (t < lifetime) ? t : lifetime;
                hit = true;
            }
        }

        if (!hit) {
            int next_cell;
            double2 intersection_point;
            bool success = GetNextCell(current_cell, 
                orbit, 
                entry_point, 
                sp->cells_per_row, 
                sp->cell_size, 
                sp->impurity_spawn_range, 
                &next_cell, 
                &intersection_point
            );

            if (!success) {
                // The End?
                break;
            }

            current_cell = next_cell;
            entry_point = intersection_point;
        }
    }

    return lifetime;
}

inline double lifetime(const double max_lifetime, const double2 pos, const double phi, BUFFER_ARGS)
{
    bool clockwise = (sp->is_clockwise == 1);

    const v2 vel = (double2)(cos(phi), sin(phi)) * sp->particle_speed;
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    const double2 center = GetCyclotronOrbitCenter(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    const Orbit orbit(center, orbit_radius, clockwise);

    // @Todo, Check if we start inside an impurity.
    int current_cell = get_cell_index(pos, sp->impurity_spawn_range, sp->max_expected_impurities_in_cell);
    int impurity_start = cell_indices[current_cell];
    // ....

    // @Todo, bound intersect?
    
    double lt = TraceOrbit(pos, phi, orbit, current_cell, sp, impurities, cell_indices);
    return lt;
}

inline double sim_phi_lifetime(const double2 pos, int quadrant, int step, BUFFER_ARGS)
{
    bool clockwise = (sp->is_clockwise == 1);
    bool incoherent = (sp->is_incoherent == 1);
    bool diag_regions = (sp->is_diag_regions == 1);

    const double phi = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;

    // Move to Orbit?
    const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
    return lifetime(min(15.0 * sp->tau, bound_time), pos, phi, sp, impurities, cell_indices);
}
