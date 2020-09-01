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
        particle_count = sp.dim * sp.dim;

        main_buffer.clear();
        main_buffer.resize(particle_count, 0);
    }

    sp = p_sp;
    first_run = false;

    return true;
}
