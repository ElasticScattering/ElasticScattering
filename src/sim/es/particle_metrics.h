#pragma once

typedef struct ParticleMetrics {
    long cells_passed;
    long impurities_tested;

    long particles_inside_impurity;
    long particles_escaped;
    long particles_at_bound;

    long max_impurities_tested;
    long max_cells_passed;

#ifndef DEVICE_PROGRAM
    ParticleMetrics() {
        cells_passed = 0;
        impurities_tested = 0;
        particles_inside_impurity = 0;
        particles_escaped = 0;
        particles_at_bound = 0;
        max_impurities_tested = 0;
        max_cells_passed = 0;
    }
#endif
} ParticleMetrics;
