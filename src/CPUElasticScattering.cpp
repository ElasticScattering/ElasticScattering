#include "ElasticScattering.h"

double CPUElasticScattering::Compute(const SimulationParameters &p_sp)
{
    PrepareCompute(p_sp);

    int limit = sp.dim - 1;

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++)
        {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            double particle_result = (sp.mode == MODE_DIR_LIFETIME) ? single_lifetime(pos, &sp, impurities) : phi_lifetime(pos, &sp, impurities);
            main_buffer[j * sp.dim + i] = particle_result * GetWeight(i, j, sp.dim);
        }
    }

    return ComputeResult(main_buffer);
}

bool CPUElasticScattering::PrepareCompute(const SimulationParameters &p_sp) {
    if (first_run || ImpuritySettingsChanged(p_sp))
        GenerateImpurities(p_sp);

    if (first_run || (sp.dim != p_sp.dim)) {
        main_buffer.clear();
        main_buffer.resize(p_sp.particle_count, 0);
    }

    sp.region_size = p_sp.region_size;
    sp.region_extends = p_sp.region_extends;
    sp.dim = p_sp.dim;
    sp.particle_speed = p_sp.particle_speed;
    sp.particle_mass = p_sp.particle_mass;
    sp.impurity_count = p_sp.impurity_count;
    sp.impurity_radius = p_sp.impurity_radius;
    sp.alpha = p_sp.alpha;
    sp.phi = p_sp.phi;
    sp.magnetic_field = p_sp.magnetic_field;
    sp.tau = p_sp.tau;
    sp.integrand_steps = p_sp.integrand_steps;
    sp.clockwise = p_sp.clockwise;

    sp.mode = p_sp.mode;
    sp.impurity_seed = p_sp.impurity_seed;

    sp.particle_count = sp.dim * sp.dim;
    sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
    sp.angular_speed = E * sp.magnetic_field / sp.particle_mass;
    sp.region_extends = sp.particle_speed * sp.tau;

    first_run = false;

    return true;
}
