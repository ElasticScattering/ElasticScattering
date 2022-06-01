/* -------------------------------------------------------------------------
    This code is part of ElasticScattering.

    Copyright(C) 2022 Elastic Scattering developers

    This program is free software : you can redistribute it and /or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#pragma once

typedef struct ParticleMetrics {
    int cells_passed;
    int impurities_tested;

    int particles_inside_impurity;
    int particles_escaped;
    int particles_at_bound;

    int max_impurities_tested;
    int max_cells_passed;

#ifndef DEVICE_PROGRAM
    ParticleMetrics() {
        cells_passed              = 0;
        impurities_tested         = 0;
        particles_inside_impurity = 0;
        particles_escaped         = 0;
        particles_at_bound        = 0;
        max_impurities_tested     = 0;
        max_cells_passed          = 0;
    }
#endif
} ParticleMetrics;
