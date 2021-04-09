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
