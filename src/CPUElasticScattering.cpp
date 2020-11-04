#include "ElasticScattering.h"
#include "src/escl/common.h"

ScatterResult CPUElasticScattering::ComputeResult(ScatteringParameters& p_sp)
{
    PrepareCompute(p_sp);

    ResultBuffer buffer(p_sp.dim * p_sp.dim);

    // GPU kernel works only with even work size.
    int limit = sp.dim - 1;

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++)
        {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            buffer.intermediate_results[j * sp.dim + i] = sim_phi_lifetime(pos, &p_sp, grid.impurities, grid.imp_index);;
        }
    }

    // Apply weights for integration.
    for (int j = 0; j < limit; j++)
        for (int i = 0; i < limit; i++) {
            buffer.intermediate_results[j * sp.dim + i].xx *= GetWeight2D(i, j, limit);
            buffer.intermediate_results[j * sp.dim + i].xy *= GetWeight2D(i, j, limit);
        }

    ScatterResult sr;
    sr = FinishResult(buffer);

    return sr;
}


bool CPUElasticScattering::ComputeSingle(ScatteringParameters& p_sp, double& result)
{
    if (!PrepareCompute(p_sp)) return false;

    std::vector<double> buffer;
    buffer.resize(p_sp.dim * p_sp.dim);
    
    // GPU kernel works only with even work size.
    int limit = sp.dim - 1;

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++)
        {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            buffer[j * sp.dim + i] = phi_lifetime(pos, &p_sp, grid.impurities, grid.imp_index);;
        }
    }

    // Apply weights for integration.
    for (int j = 0; j < limit; j++)
        for (int i = 0; i < limit; i++)
            buffer[j * sp.dim + i] *= GetWeight2D(i, j, limit);
    
    result = FinishSingle(buffer);

    return true;
}

bool CPUElasticScattering::PrepareCompute(ScatteringParameters &p_sp) {
    CompleteSimulationParameters(p_sp);
    
    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed  = (sp.dim != p_sp.dim);

    if (!first_run && !impurities_changed && !work_size_changed) return false;

    sp = p_sp;
        
    if (first_run || impurities_changed)
        grid.Generate(sp);

    first_run = false;

    return true;
}
