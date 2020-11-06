#pragma once

#include "device_macros.h"
#include "details.h"
#include "impurity_grid.h"

#include "windows.h"

double SIM_lifetime(const v2 pos, const double phi, ScatteringParameters* sp, std::vector<v2>& impurities, std::vector<int>& cell_index) {
    const v2 unit = { cos(phi), sin(phi) };
    const v2 vel = unit * sp->particle_speed;

    double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    v2 range = { -sp->region_extends, sp->region_size + sp->region_extends };

    int2 cell = to_grid(pos.x, pos.y, range, 113);

    int index = cell.y * sp->max_expected_impurities_in_cell + cell.x;
    int starting_cell_index = cell_index[index];
    int ending_cell_index = cell_index[index + 1];

    double max_lifetime = 15.0 * sp->tau;
    double lifetime_old = max_lifetime;
    bool hit = false;
    while (!hit)
    {
        for (int i = starting_cell_index; i < ending_cell_index; i++)
        {
            double lt = 0; // lifetime_old(max_life(pos, impurities[i], impurity_radius_sq, unit, vel);

            lifetime_old = (lt < lifetime_old) ? lt : lifetime_old;

            if (lifetime_old == 0) {
                break;
            }
        }

        hit = lifetime_old < max_lifetime;
        if (!hit) {
            // calculate next grid cell to move to...
            // index = ....
            // starting_index / ending_index
        }
    }

    return lifetime_old;
}
inline double lifetime(const double max_lifetime, const double2 pos, const double phi, const bool clockwise, BUFFER_ARGS)
{
    const v2 unit = { cos(phi), sin(phi) };
    const v2 vel = unit * sp->particle_speed;
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    const double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);

    v2 range = { -sp->region_extends, sp->region_size + sp->region_extends };

    int2 cell = to_grid(pos.x, pos.y, range, 113);

    int index = cell.y * sp->max_expected_impurities_in_cell + cell.x;
    int starting_cell_index = cell_index[index];
    int ending_cell_index = cell_index[index + 1];
    
    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    //max_lifetime = 15.0 * sp->tau;  Move to SP
    double lifetime_old = max_lifetime;
    bool hit = false;

    //...

    return lifetime_old;
}

inline double trace_grid()
{
    /*
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
    */
}

// ! Deprecated
inline double lifetime_old(const double max_lifetime, const double2 pos, const double phi, const bool clockwise, BUFFER_ARGS)
{
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = { cos(phi), sin(phi) };
    vel = vel * sp->particle_speed;
    const double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);
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