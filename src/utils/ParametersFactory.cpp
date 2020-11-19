#include "ParametersFactory.h"
#include "src/scattering/escl/constants.h"

ScatteringParameters& ParametersFactory::GenerateSimulation() {
    ScatteringParameters sp;

    sp.integrand_steps = 13;
    sp.dim = 128;

    sp.temperature = 4;
    sp.tau = 1e-11;
    sp.magnetic_field = 0;
    sp.alpha = 0.3;
    sp.particle_speed = 1.67834e5;
    
    sp.impurity_density = 5.34e14;
    sp.impurity_radius = 1.11e-8;
    sp.impurity_count = 100;
    sp.impurity_seed = 0;
    sp.region_extends = 1e-6;
    sp.region_size = 4e-6;

    sp.is_diag_regions = 0;
    sp.is_clockwise = 0;
    sp.is_incoherent = 1;

    return sp;
}

ScatteringParameters& ParametersFactory::GenerateDefault() {
	ScatteringParameters sp;

    sp.integrand_steps = 13;
    sp.dim = 128;

    sp.temperature = 4;
    sp.tau = 1.5e-12;
    sp.magnetic_field = 0;
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

    return sp;
}

ScatteringParameters& ParametersFactory::GenerateNoImpurities() {
    ScatteringParameters sp;

    sp.integrand_steps = 13;
    sp.dim = 128;
    
    sp.temperature = 4;
    sp.tau = 1.5e-12;
    sp.magnetic_field = 0;
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

    return sp;
}

ScatteringParameters& ParametersFactory::GenerateMinimal() {
    ScatteringParameters sp;

    sp.integrand_steps = 7;
    sp.dim = 64;

    sp.temperature = 4;
    sp.tau = 1.5e-12;
    sp.magnetic_field = 0;
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


    return sp;
}
