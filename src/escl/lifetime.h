#pragma once

#include "device_macros.h"
#include "details.h"
#include "impurity_grid.h"

#include "windows.h"

//Backwards compatible
inline double lifetimeB(const double max_lifetime, const double2 pos, const double phi, const bool clockwise, BUFFER_ARGS)
{
    const double orbit_radius = sp->particle_speed / sp->angular_speed;
    double2 vel = { cos(phi), sin(phi) };
    vel = vel * sp->particle_speed;
    const double2 center = GetCyclotronOrbit(pos, vel, orbit_radius, sp->particle_speed, clockwise);
    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = max_lifetime;

    for (int i = 0; i < sp->impurity_count; i++) {
        const double2 impurity = impurities[i];

        double2 d = pos - impurity;
        if (impurity_radius_sq > dot(d, d))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, orbit_radius, impurity, sp->impurity_radius))
        {
            double t = GetFirstCrossTime(center, pos, impurity, orbit_radius, sp->impurity_radius, sp->angular_speed, clockwise);

            lifetime = (t < lifetime) ? t : lifetime;
        }
    }

    return lifetime;
}

//Backwards compatible
inline double lifetime0(const double2 pos, const double phi, BUFFER_ARGS)
{
    const double2 unit = { cos(phi), sin(phi) };
    const double2 vel = unit * sp->particle_speed;

    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double lifetime = 15.0 * sp->tau;

    for (int i = 0; i < sp->impurity_count; i++) {
        const double2 imp_pos = impurities[i];

        double lt = ComputeLifetime0(pos, imp_pos, impurity_radius_sq, unit, vel);

        lifetime = (lt < lifetime) ? lt : lifetime;
        if (lifetime == 0) {
            break;
        }
    }

    return lifetime;
}

inline double ComputeLifetime0(double2 pos, double2 imp_pos, double impurity_radius_sq, double2 unit, double2 vel)
{
    const double inner = (imp_pos.x - pos.x) * unit.x + (imp_pos.y - pos.y) * unit.y;
    const double2 projected = pos + unit * inner;

    const double2 d = projected - imp_pos;
    const double diff = impurity_radius_sq - dot(d, d);
    if (diff < 0.0) {
        return DBL_MAX;
    }

    double L = sqrt(diff);

    double2 time_taken;
    double min_time_taken;

    if (fabs(vel.x) > (fabs(vel.y) * 1e-9)) {
        time_taken.x = -((projected.x - L * unit.x) - pos.x) / vel.x;
        time_taken.y = -((projected.x + L * unit.x) - pos.x) / vel.x;
    }
    else {
        time_taken.x = -((projected.y - L * unit.y) - pos.y) / vel.y;
        time_taken.y = -((projected.y + L * unit.y) - pos.y) / vel.y;
    }

    if ((time_taken.x * time_taken.y) < 0) {
        return 0;
    }

    if (time_taken.y > 0 && time_taken.x > 0) {
        return min(time_taken.x, time_taken.y);
    }
    else if (time_taken.y > 0) {
        return time_taken.y;
    }
    else {
        return time_taken.x;
    }
}

double lifetime0C(const double2 pos, const double phi, BUFFER_ARGS)
{
    const double2 unit = { cos(phi), sin(phi) };
    const double2 vel = unit * sp->particle_speed;

    const double impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;

    double2 range = { -sp->region_extends, sp->region_size + sp->region_extends };

    int2 cell = to_grid(pos.x, pos.y, range, 113);

    int index = cell.y * sp->max_expected_impurities_in_cell + cell.x;
    int starting_cell_index = cell_index[index];
    int ending_cell_index = cell_index[index + 1];

    double max_lifetime = 15.0 * sp->tau;
    double lifetime = max_lifetime;
    bool hit = false;
    while (!hit)
    {
        for (int i = starting_cell_index; i < ending_cell_index; i++)
        {
            double lt = ComputeLifetime0(pos, impurities[i], impurity_radius_sq, unit, vel);

            lifetime = (lt < lifetime) ? lt : lifetime;
            
            if (lifetime == 0) {
                break;
            }
        }

        hit = lifetime < max_lifetime;
        if (!hit) {
            // calculate next grid cell to move to...
            // index = ....
            // starting_index / ending_index
        }
    }

    return lifetime;
}