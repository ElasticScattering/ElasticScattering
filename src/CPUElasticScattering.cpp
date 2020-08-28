#include <windows.h>

#include <random>
#include <limits>
#include <assert.h>
#include <iostream>

#include "ElasticScattering.h"
#include "Details.h"

void CPUElasticScattering::Init(bool show_info)
{
}

void CPUElasticScattering::PrepareCompute(const SimulationParameters *p_sp) {
    sp = new SimulationParameters;
    sp->region_size = p_sp->region_size;
    sp->dim = p_sp->dim;
    sp->particle_speed = p_sp->particle_speed;
    sp->particle_mass = p_sp->particle_mass;
    sp->impurity_count = p_sp->impurity_count;
    sp->impurity_radius = p_sp->impurity_radius;
    sp->alpha = p_sp->alpha;
    sp->phi = p_sp->phi;
    sp->magnetic_field = p_sp->magnetic_field;
    sp->tau = p_sp->tau;
    sp->integrand_steps = p_sp->integrand_steps;
    sp->clockwise = p_sp->clockwise;

    sp->particle_count = sp->dim * sp->dim;
    sp->impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;
    sp->angular_speed = E * sp->magnetic_field / sp->particle_mass;

    if (sp->angular_speed == 0) {
        sp->particle_max_lifetime = sp->tau;
    }
    else {
        double bound_time = GetBoundTime(sp->phi, sp->alpha, sp->angular_speed, false, false); // @todo, electorn/clockwise from parameters.
        sp->particle_max_lifetime = MIN(sp->tau, bound_time);
    }

    if (false) {
        std::cout << "Particle max lifetime: " << sp->particle_max_lifetime << std::endl;
    }

    // Initialize arrays.
    impurities.clear();
    impurities.resize(sp->impurity_count);

    if (false)
        std::cout << "Impurity region: " << -sp->particle_speed * sp->tau << ", " << sp->region_size + sp->particle_speed * sp->tau << std::endl;

    std::uniform_real_distribution<double> unif(-sp->particle_speed * sp->tau, sp->region_size + sp->particle_speed * sp->tau);
    std::random_device r;
    std::default_random_engine re(0);

    for (int i = 0; i < sp->impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    main_buffer.clear();
    main_buffer.resize(sp->particle_count, 0);
}

double CPUElasticScattering::Compute(Mode p_mode, const SimulationParameters* p_sp)
{
    PrepareCompute(p_sp);
    mode = p_mode;

    double total_time;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    //std::cout << "Simulating elastic scattering on the CPU..." << std::endl;

    int limit = sp->dim - 1;
    QueryPerformanceCounter(&beginClock);
    if (mode == Mode::AVG_LIFETIME)
    {
        for (int j = 0; j < limit; j++)
            for (int i = 0; i < limit; i++)
            {
                v2 pos;
                pos.x = sp->region_size * (double(i) / (sp->dim - 2));
                pos.y = sp->region_size * (double(j) / (sp->dim - 2));

                double lifetime = sp->particle_max_lifetime;

                v2 vel = { sp->particle_speed * cos(sp->phi), sp->particle_speed * sin(sp->phi) };

                double res = (sp->angular_speed == 0) ? ComputeA(pos, vel) : ComputeB(pos, vel);

                bool is_edge = (i == 0) || (i == sp->dim - 2) || (j == 0) || (j == sp->dim - 2);

                double w = 1.0;
                if (!is_edge)
                {
                    w  = ((i % 2) == 0) ? 2.0 : 4.0;
                    w *= ((j % 2) == 0) ? 2.0 : 4.0;
                }

                main_buffer[j * sp->dim + i] = w * res;
            }
    }
    else if (mode == Mode::SIGMA_XX)
    {
        for (int y = 0; y < limit; y++)
        {
            for (int x = 0; x < limit; x++)
            {
                v2 pos;
                pos.x = sp->region_size * (double(x) / (sp->dim - 2));
                pos.y = sp->region_size * (double(y) / (sp->dim - 2));

                double lifetime = sp->particle_max_lifetime;

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

                        v2 vel = { sp->particle_speed * cos(sp->phi), sp->particle_speed * sin(sp->phi) };

                        double lt = (sp->angular_speed == 0) ? ComputeA(pos, vel) : ComputeB(pos, vel);

                        double z = exp(-lt / sp->tau);

                        double r = cos(sp->phi) - cos(sp->phi) * z;
                        r += sp->angular_speed * sp->tau * sin(sp->phi) * z;
                        r -= sp->angular_speed * sp->tau * sin(sp->phi);
                        r *= sp->tau;

                        double rxx = r * cos(sp->phi);
                        double rxy = r * sin(sp->phi);

                        bool edge_item = (i == 0 || i == sp->dim - 1);

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
                    w  = ((x % 2) == 0) ? 2.0 : 4.0;
                    w *= ((y % 2) == 0) ? 2.0 : 4.0;
                }

                main_buffer[y * sp->dim + x] = w * integral;
            }
        }
    }

    double total = 0;
    for (int j = 0; j < limit; j++)
        for (int i = 0; i < limit; i++)
            total += main_buffer[j*sp->dim + i];

    int actual_particle_count = (sp->dim - 1) * (sp->dim - 1);

    QueryPerformanceCounter(&endClock);

    /*
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;
    
    MakeTexture(sp);
    */
    
    return total / (double)actual_particle_count;
}

double CPUElasticScattering::ComputeA(const v2 pos, const v2 vel)
{
    const v2 unit = { cos(sp->phi), sin(sp->phi) };
    double lifetime = sp->tau;

    for (int k = 0; k < sp->impurity_count; k++)
    {
        const v2 ip = impurities[k];
        const double inner = (ip.x - pos.x) * unit.x + (ip.y - pos.y) * unit.y;
        const v2 projected = { pos.x + inner * unit.x, pos.y + inner * unit.y };

        const double diffx = projected.x - ip.x;
        const double diffy = projected.y - ip.y;
        const double a = diffx*diffx + diffy*diffy;
        if (a > sp->impurity_radius_sq) {
            continue; //@Speedup, if distance is greater than current min continue as well?
        }

        double L = sqrt(sp->impurity_radius_sq - a);

        v2 time_taken;
        if (vel.x != 0) time_taken = { -((projected.x - L * unit.x) - pos.x) / vel.x, -((projected.x + L * unit.x) - pos.x) / vel.x };
        else            time_taken = { -((projected.y - L * unit.y) - pos.y) / vel.y, -((projected.y + L * unit.y) - pos.y) / vel.y };

        if ((time_taken.x * time_taken.y) < 0)
            return 0;

        if (time_taken.x > 0 && time_taken.x < lifetime) {
            lifetime = time_taken.x;
        }
        if (time_taken.y > 0 && time_taken.y < lifetime) {
            lifetime = time_taken.y;
        }
    }

    return lifetime;
}

double CPUElasticScattering::ComputeB(const v2 pos, const v2 vel)
{
    double vf = sp->particle_speed; // sqrt(vel.x * vel.x + vel.y * vel.y);
    double radius = vf / sp->angular_speed;
    auto center = GetCyclotronOrbit(pos, vel, radius, vf, sp->clockwise);

    double lifetime = sp->particle_max_lifetime;
    for (int k = 0; k < sp->impurity_count; k++)
    {
        const v2 ip = impurities[k];

        if (sp->impurity_radius_sq > pow(pos.x - ip.x, 2) + pow(pos.y - ip.y, 2))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, radius, ip, sp->impurity_radius))
        {
            double t = GetFirstCrossTime(center, pos, ip, radius, sp->impurity_radius, sp->angular_speed, sp->clockwise);

            assert(t >= 0);

            if (t < lifetime)
                lifetime = t;
        }
    }

    return lifetime;
}
