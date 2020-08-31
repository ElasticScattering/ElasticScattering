#include <windows.h>

#include <assert.h>
#include <iostream>

#include "ElasticScattering.h"

void CPUElasticScattering::PrepareCompute(const SimulationParameters *p_sp) {
    bool first_run = (last_sp == nullptr);

    bool need_impurity_update = first_run || (sp->impurity_count != last_sp->impurity_count || sp->region_extends != last_sp->region_extends || sp->region_size != last_sp->region_size);
    bool need_particle_update = first_run || (sp->dim != last_sp->dim);

    sp = new SimulationParameters;
    sp->region_size        = p_sp->region_size;
    sp->region_extends     = p_sp->region_extends;
    sp->dim                = p_sp->dim;
    sp->particle_speed     = p_sp->particle_speed;
    sp->particle_mass      = p_sp->particle_mass;
    sp->impurity_count     = p_sp->impurity_count;
    sp->impurity_radius    = p_sp->impurity_radius;
    sp->alpha              = p_sp->alpha;
    sp->phi                = p_sp->phi;
    sp->magnetic_field     = p_sp->magnetic_field;
    sp->tau                = p_sp->tau;
    sp->integrand_steps    = p_sp->integrand_steps;
    sp->clockwise          = p_sp->clockwise;
    
    sp->mode               = p_sp->mode;
    sp->impurity_seed      = p_sp->impurity_seed;

    sp->particle_count     = sp->dim * sp->dim;
    sp->impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;
    sp->angular_speed      = E * sp->magnetic_field / sp->particle_mass;
    sp->region_extends     = sp->particle_speed * sp->tau;

    if (need_impurity_update)
        GenerateImpurities();

    if (need_particle_update) {
        main_buffer.clear();
        main_buffer.resize(sp->particle_count, 0);
    }
}

double CPUElasticScattering::Compute(const SimulationParameters* p_sp)
{
    PrepareCompute(p_sp);

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    QueryPerformanceCounter(&beginClock);
    
    int limit = sp->dim - 1;
    auto imps = impurities.data();

    for (int j = 0; j < limit; j++)
        for (int i = 0; i < limit; i++)
        {
            v2 pos(i, j);
            pos = pos * (sp->region_size / (double)(sp->dim - 2));

            double particle_result = (sp->mode == MODE_DIR_LIFETIME) ? single_lifetime(pos, sp, imps) : phi_lifetime(pos, sp, imps);
            main_buffer[j * sp->dim + i] = particle_result * GetWeight(i, j, sp->dim);
        }

    double result = ComputeResult(main_buffer);

    QueryPerformanceCounter(&endClock);

    //double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    //std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;
    
    return result;
}
