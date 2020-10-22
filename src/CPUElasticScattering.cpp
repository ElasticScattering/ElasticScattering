#include "ElasticScattering.h"
#include "src/escl/common.h"

bool CPUElasticScattering::Compute(ScatteringParameters& p_sp, double& result)
{
    if (!PrepareCompute(p_sp)) return false;

    // GPU kernel works only with even work size.
    int limit = sp.dim - 1;

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++)
        {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            main_buffer[j * sp.dim + i] = (sp.mode == MODE_DIR_LIFETIME) ? single_lifetime(pos, sp.phi, &sp, grid.impurities, grid.imp_index) : 
                                                                           phi_lifetime(pos, &sp, grid.impurities, grid.imp_index);
        }
    }

    // Apply weights for integration.
    for (int j = 0; j < limit; j++)
        for (int i = 0; i < limit; i++)
            main_buffer[j * sp.dim + i] *= GetWeight2D(i, j, limit);
    
    result = ComputeResult(main_buffer);

    return true;
}

bool CPUElasticScattering::PrepareCompute(ScatteringParameters &p_sp) {
    CompleteSimulationParameters(p_sp);
    
    if (!first_run && (p_sp == sp)) return false;

    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed  = (sp.dim != p_sp.dim);

    sp = p_sp;
        
    if (first_run || impurities_changed)
        grid.Generate(sp);

    if (first_run || work_size_changed) {
        main_buffer.clear();
        main_buffer.resize(particle_count, 0);
    }

    first_run = false;

    return true;
}
