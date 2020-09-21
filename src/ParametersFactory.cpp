#include "ParametersFactory.h"
#include "src/escl/constants.h"

SimulationParameters ParametersFactory::GenerateDefault() {
	SimulationParameters sp;

    sp.integrand_steps = 13;
    sp.dim = 128;

    sp.temperature = 4;
    sp.tau = 1.5e-12;
    sp.magnetic_field = 0;
    sp.phi = 1.0;
    sp.alpha = PI / 4.0;
    sp.particle_speed = 1.67834e5;
    sp.impurity_density = 1e12;

    sp.impurity_count = 100;
    sp.impurity_radius = 2e-9;
    sp.impurity_seed = 0;
    sp.region_size = 1e-5;
    sp.region_extends = sp.particle_speed * sp.tau * 15.0;

    sp.is_diag_regions = 0;
    sp.is_clockwise = 0;
    sp.is_incoherent = 1;

    sp.mode = MODE_DIR_LIFETIME;

    return sp;
}

SimulationParameters ParametersFactory::GenerateNoImpurities() {
    SimulationParameters sp;

    sp.integrand_steps = 13;
    sp.dim = 128;
    
    sp.temperature = 4;
    sp.tau = 1.5e-12;
    sp.magnetic_field = 0;
    sp.phi = 1.0;
    sp.alpha = PI / 4.0;
    sp.particle_speed = 1.67834e5;

    sp.impurity_count = 1;
    sp.impurity_radius = 1e-16;
    sp.impurity_seed = 0;
    sp.region_size = 1e-5;
    sp.region_extends = sp.particle_speed * sp.tau * 15.0;
    sp.impurity_density = 1e12;

    sp.is_diag_regions = 0;
    sp.is_clockwise = 0;
    sp.is_incoherent = 1;

    sp.mode = MODE_DIR_LIFETIME;

    return sp;
}

SimulationParameters ParametersFactory::GenerateMinimal() {
    SimulationParameters sp;

    sp.integrand_steps = 7;
    sp.dim = 64;

    sp.temperature = 4;
    sp.tau = 1.5e-12;
    sp.magnetic_field = 0;
    sp.phi = 1.0;
    sp.alpha = PI / 4.0;
    sp.particle_speed = 1.67834e5;

    sp.impurity_count = 1;
    sp.impurity_radius = 2e-9;
    sp.impurity_seed = 0;
    sp.region_size = 1e-5;
    sp.region_extends = sp.particle_speed * sp.tau * 15.0;
    sp.impurity_density = 1e1;

    sp.is_diag_regions = 0;
    sp.is_clockwise = 0;
    sp.is_incoherent = 1;

    sp.mode = MODE_DIR_LIFETIME;

    return sp;
}
