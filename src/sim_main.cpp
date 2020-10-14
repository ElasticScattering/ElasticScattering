#pragma once

#include "sim_main.h"

#include "escl/constants.h"
#include "utils/ParametersFactory.h"
#include "datastructures/SimulationConfiguration.h"
#include "datastructures/SimulationResult.h"
#include "Logger.h"
#include <windows.h>

int sim_main(const InitParameters& init)
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    ScatteringParameters sp = ParametersFactory::GenerateSimulation();

    SimulationConfiguration sim_params;
    sim_params.runs = 10;
    sim_params.samples_per_run = 2;
    sim_params.magnetic_field_min = 0.01;
    sim_params.magnetic_field_max = 40;
    sim_params.scattering_params = sp;

    const std::vector<double> temperatures{ 15, 60 };
    PrintInfo(sim_params, temperatures.size());

    SimulationElasticScattering es(init);

    for (int i = 0; i < temperatures.size(); i++) {
        sp.temperature = temperatures[i];

        SimulationResult sr(sim_params.runs);

        QueryPerformanceCounter(&beginClock);
        ComputeIteration(es, sim_params, sr);
        QueryPerformanceCounter(&endClock);
        sr.time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        Logger::LogResult(sim_params, sr);
        printf("Simulation completed (%d/%d)\n", i + 1, temperatures.size());
    }

    return 0;
}

void ComputeIteration(SimulationElasticScattering &es, SimulationConfiguration& sp, SimulationResult& sr)
{
    double coherent_tau = sp.scattering_params.tau;
    bool run_incoherent = sp.scattering_params.alpha > 0.000001;
    bool run_coherent = abs(sp.scattering_params.alpha - (PI / 4)) > 0.000001;
    
    double step_size = (sp.magnetic_field_max - sp.magnetic_field_min) / sp.runs;

    printf("magnetic_field sigma_xx_inc sigma_xx_coh sigma_xy_inc sigma_xy_coh delta_xx\n");

    for (int i = 0; i < sp.runs; i++) {
        sp.scattering_params.magnetic_field = sp.magnetic_field_min + step_size * i;

        {
            double total_sigma_xx = 0;
            double total_sigma_xxi = 0;
            double total_sigma_xx_sq = 0;

            double total_sigma_xy = 0;
            double total_sigma_xyi = 0;
            double total_sigma_xy_sq = 0;

            for (int j = 0; j < sp.samples_per_run; j++) {
                sp.scattering_params.impurity_seed = 1123 + j * 831;
                v2 result, resulti;

                if (run_incoherent) {
                    sp.scattering_params.is_incoherent = 1;
                    es.Compute(sp.scattering_params, resulti);
                }
                else {
                    resulti = 0;
                }

                if (run_coherent) {
                    sp.scattering_params.is_incoherent = 0;
                    sp.scattering_params.tau = coherent_tau;
                    es.Compute(sp.scattering_params, result);
                }
                else {
                    result = 0;
                }

                double sxx = result.x / 1e8;
                double sxy = result.y / 1e8;
                
                double sxx_i = resulti.x / 1e8;
                double sxy_i = resulti.y / 1e8;

                total_sigma_xx += sxx;
                total_sigma_xxi += sxx_i;
                total_sigma_xx_sq += (sxx_i + sxx) * (sxx_i + sxx);

                total_sigma_xy += sxy;
                total_sigma_xyi += sxy_i;
            }

            Result r;
            r.x = sp.scattering_params.magnetic_field;
            //sr.xs_temperature[i] = sp.scattering_params.temperature;

            {
                double sxx_sq_exp = total_sigma_xx_sq / (double)(sp.samples_per_run) + 1e-15;
                double sxx_exp = (total_sigma_xx + total_sigma_xxi) / (double)(sp.samples_per_run);

                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(sp.samples_per_run - 1));
                
                r.xx  = total_sigma_xx / (double)(sp.samples_per_run);
                r.xxi = total_sigma_xxi / (double)(sp.samples_per_run);

                r.xxd = sxx_std / sxx_exp;

                r.xy  = total_sigma_xy / (double)sp.samples_per_run;
                r.xyi = total_sigma_xyi / (double)sp.samples_per_run;
            }

            sr.results[i] = r;
            printf("%f %e %e %e %e %f\n", r.x, r.xxi, r.xx, r.xyi, r.xy, r.xxd);
        }
    }
}

void PrintInfo(const SimulationConfiguration& sp, int count)
{
    const int imp_count = sp.scattering_params.impurity_density * pow(sp.scattering_params.region_extends + sp.scattering_params.region_size, 2);
    const long long intersects_each = sp.samples_per_run * imp_count * 2.0 * pow(sp.scattering_params.dim, 2) * sp.scattering_params.integrand_steps * 4.0;
    const long long intersects = sp.runs * intersects_each;
    printf("Simulation info:\n");
    printf("Total intersections: %e\n", intersects * (long long)count);
    printf("Time estimate per data point: %f minutes.\n", (float)intersects_each / 2e9 / 60);
    printf("Time estimate per temperature: %f hours \n", (float)intersects / 2e9 / 3600.0);
    printf("Time estimate total: %f hours \n\n", (float)intersects / 2e9 / 3600.0 * count);
}
