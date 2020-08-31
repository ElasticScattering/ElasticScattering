#include <windows.h>

#include <assert.h>
#include <iostream>

#include "ElasticScattering.h"

void CPUElasticScattering::Init(bool show_info)
{
}

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

void CPUElasticScattering::SigmaXX()
{
    int limit = sp->dim - 1;
    auto imps = impurities.data();

    for (int y = 0; y < limit; y++)
    {
        for (int x = 0; x < limit; x++)
        {
            v2 pos;
            pos.x = sp->region_size * (x / (double)(sp->dim - 2));
            pos.y = sp->region_size * (y / (double)(sp->dim - 2));

            double angle_area = sp->alpha * 2.0;
            double step_size = angle_area / (sp->integrand_steps - 1);
            double integral = 0;

            for (int j = 0; j < 4; j++)
            {
                double start = -sp->alpha + j * (PI * 0.5);
                double total = 0.0;

                for (int i = 0; i < sp->integrand_steps; i++)
                {
                    sp->phi = start + i * step_size;

                    double particle_lifetime;
                    if (sp->angular_speed != 0) {
                        bool clockwise = (sp->clockwise == 1);
                        if (clockwise) sp->angular_speed *= -1;

                        double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, clockwise, false);
                        particle_lifetime = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, imps);
                    }
                    else {
                        particle_lifetime = lifetime0(pos, sp, imps);
                    }

                    double z = exp(-particle_lifetime / sp->tau);

                    double r = cos(sp->phi) - cos(sp->phi + sp->angular_speed * particle_lifetime) * z;
                    r += sp->angular_speed * sp->tau * sin(sp->phi + sp->angular_speed * particle_lifetime) * z;
                    r -= sp->angular_speed * sp->tau * sin(sp->phi);
                    r *= sp->tau;

                    double rxx = r * cos(sp->phi);
                    double rxy = r * sin(sp->phi);

                    bool edge_item = (i == 0 || i == sp->integrand_steps - 1);

                    double w = 1.0;

                    if (!edge_item) {
                        w = ((i % 2) == 0) ? 2.0 : 4.0;
                    }

                    total += rxx * w;
                }

                integral += total * angle_area / ((sp->integrand_steps - 1) * 3.0);
            }

            bool is_edge = (x == 0) || (x == (sp->dim - 2)) || (y == 0) || (y == (sp->dim - 2));

            double w = 1.0;
            if (!is_edge)
            {
                w = ((x % 2) == 0) ? 2.0 : 4.0;
                w *= ((y % 2) == 0) ? 2.0 : 4.0;
            }

            main_buffer[y * sp->dim + x] = w * integral;
        }
    }
}

void CPUElasticScattering::Lifetime()
{
    int limit = sp->dim - 1;

    auto imps = impurities.data();

    for (int j = 0; j < limit; j++)
        for (int i = 0; i < limit; i++)
        {
            v2 pos;
            pos.x = sp->region_size * (i / (double)(sp->dim - 2));
            pos.y = sp->region_size * (j / (double)(sp->dim - 2));

            double particle_lifetime;
            if (sp->angular_speed != 0) {
                bool clockwise = (sp->clockwise == 1);

                double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, clockwise, false);
                particle_lifetime = lifetimeB(min(sp->tau, bound_time), pos, clockwise, sp, imps);
            }
            else {
                particle_lifetime = lifetime0(pos, sp, imps);
            }

            bool is_edge = (i == 0) || (i == sp->dim - 2) || (j == 0) || (j == sp->dim - 2);

            double w = 1.0;
            if (!is_edge)
            {
                w = ((i % 2) == 0) ? 2.0 : 4.0;
                w *= ((j % 2) == 0) ? 2.0 : 4.0;
            }

            main_buffer[j * sp->dim + i] = w * particle_lifetime;
        }
}

double CPUElasticScattering::Compute(Mode p_mode, const SimulationParameters* p_sp)
{
    PrepareCompute(p_sp);
    mode = p_mode;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    QueryPerformanceCounter(&beginClock);
    
    if      (mode == Mode::AVG_LIFETIME) Lifetime();
    else if (mode == Mode::SIGMA_XX)     SigmaXX();

    double result = ComputeResult(main_buffer);

    QueryPerformanceCounter(&endClock);

    //double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    //std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;
    
    return result;
}