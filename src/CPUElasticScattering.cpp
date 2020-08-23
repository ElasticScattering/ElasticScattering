#include <windows.h>

#include <random>
#include <limits>
#include <assert.h>
#include <iostream>

#include "ElasticScattering.h"
#include "Details.h"

void CPUElasticScattering::Init(InitParameters p_ip, SimulationParameters p_sp)
{
    sp = p_sp;
    mode = p_ip.mode;

    if (sp.angular_speed == 0) sp.particle_max_lifetime = sp.tau;
    else {
        double bound_time = GetBoundTime(sp.phi, sp.alpha, sp.angular_speed, true, false);
        sp.particle_max_lifetime = MIN(sp.tau, bound_time);
    }

    if (p_ip.show_info) {
        std::cout << "Particle max lifetime: " << sp.particle_max_lifetime << std::endl;
    }

    // Initialize arrays.
    impurities.clear();
    impurities.resize(sp.impurity_count);

    if (p_ip.show_info)
        std::cout << "Impurity region: " << -sp.particle_speed * sp.tau << ", " << sp.region_size + sp.particle_speed * sp.tau << std::endl;

    std::uniform_real_distribution<double> unif(-sp.particle_speed * sp.tau, sp.region_size + sp.particle_speed * sp.tau);
    std::random_device r;
    std::default_random_engine re(0);

    for (int i = 0; i < sp.impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    main_buffer.clear();
    main_buffer.resize(sp.particle_count, 0);
}

double CPUElasticScattering::Compute()
{
    double total_time;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    //std::cout << "Simulating elastic scattering on the CPU..." << std::endl;

    QueryPerformanceCounter(&beginClock);
    if (mode == Mode::AVG_LIFETIME)
    {
        for (int j = 0; j < sp.dim; j++)
            for (int i = 0; i < sp.dim; i++)
            {
                v2 pos;
                pos.x = sp.region_size * (double(i) / sp.dim);
                pos.y = sp.region_size * (double(j) / sp.dim);

                double lifetime = sp.particle_max_lifetime;

                v2 vel = { sp.particle_speed * cos(sp.phi), sp.particle_speed * sin(sp.phi) };

                double res = (sp.angular_speed == 0) ? ComputeA(pos, vel, sp) : ComputeB(pos, vel, sp);

                bool is_padding = (i == sp.dim - 1) || (j == sp.dim - 1);
                bool is_edge = (i == 0) || (i == sp.dim - 2) || (j == 0) || (j == sp.dim - 2);

                double w = is_padding ? 0.0 : 1.0;
                if (!is_edge)
                {
                    w  = ((i % 2) == 0) ? 2.0 : 4.0;
                    w *= ((j % 2) == 0) ? 2.0 : 4.0;
                }

                main_buffer[j * sp.dim + i] = w * res;
            }

        double total = 0;
        for (int i = 0; i < sp.particle_count; i++) {
            total += main_buffer[i];
        }
        
        int actual_particle_count = (sp.dim - 1) * (sp.dim - 1);
        return total / (double)actual_particle_count;
    }
    
    QueryPerformanceCounter(&endClock);
    
    return -99999;

    /*
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;
    
    MakeTexture(sp);
    */
}

double CPUElasticScattering::ComputeA(const v2 pos, const v2 vel, const SimulationParameters sp)
{
    const v2 unit = { cos(sp.phi), sin(sp.phi) };
    double lifetime = sp.tau;

    for (int k = 0; k < sp.impurity_count; k++)
    {
        const v2 ip = impurities[k];
        const double inner = (ip.x - pos.x) * unit.x + (ip.y - pos.y) * unit.y;
        const v2 projected = { pos.x + inner * unit.x, pos.y + inner * unit.y };

        const double diffx = projected.x - ip.x;
        const double diffy = projected.y - ip.y;
        const double a = diffx*diffx + diffy*diffy;
        if (a > sp.impurity_radius_sq) {
            continue; //@Speedup, if distance is greater than current min continue as well?
        }

        double L = sqrt(sp.impurity_radius_sq - a);

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

double CPUElasticScattering::ComputeB(const v2 pos, const v2 vel, const SimulationParameters sp)
{
    const bool clockwise = true;

    double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
    double radius = vf / sp.angular_speed;
    auto center = GetCyclotronOrbit(pos, vel, radius, vf, clockwise);

    double lifetime = sp.particle_max_lifetime;
    for (int k = 0; k < sp.impurity_count; k++)
    {
        const v2 ip = impurities[k];

        if (sp.impurity_radius_sq > pow(pos.x - ip.x, 2) + pow(pos.y - ip.y, 2))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, radius, ip, sp.impurity_radius))
        {
            double t = GetFirstCrossTime(center, pos, ip, radius, sp.impurity_radius, sp.angular_speed, clockwise);

            assert(t >= 0);

            if (t < lifetime)
                lifetime = t;
        }
    }

    return lifetime;
}
