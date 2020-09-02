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
            main_buffer[j * sp.dim + i] = particle_result * GetWeight2D(i, j, sp.dim);
        }
    }

    return ComputeResult(main_buffer);
}

bool CPUElasticScattering::PrepareCompute(const SimulationParameters &p_sp) {
    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed  = (sp.dim != p_sp.dim);

    sp = p_sp;
    CompleteSimulationParameters();
        
    if (first_run || impurities_changed)
        GenerateImpurities(sp);

    if (first_run || work_size_changed) {
        main_buffer.clear();
        main_buffer.resize(particle_count, 0);
    }

    first_run = false;

    return true;
}
