#pragma once
#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"
#include "utils/ErrorMacros.h"
#include "src/ParametersFactory.h"
#include "src/Logger.h"
#include "src/escl/constants.h"
#include "src/SimulationResult.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include "math.h"
#include <thread>
#include <chrono>

LARGE_INTEGER beginClock, endClock, clockFrequency;

SimulationParameters sp;

Logger logger;

void ComputeSimulation(ElasticScattering& es, SimulationResult& sr)
{
    /*
    std::vector<double> zs{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                    1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    sr.n_runs = zs.size();
    */

    sr.n_runs = 20;
    sr.xs.resize(sr.n_runs);
    sr.xs_temperature.resize(sr.n_runs);

    sr.results_xx.resize(sr.n_runs);
    sr.results_xxi.resize(sr.n_runs);
    sr.delta_xxi.resize(sr.n_runs);
    sr.results_xy.resize(sr.n_runs);
    sr.results_xyi.resize(sr.n_runs);

    sr.iterations_per_run = 3;
    sr.x_is_temperature = false;

    v2 range = { 0.01, 40 };
    double step_size = (range.y - range.x) / sr.n_runs;

    double coherent_tau = sp.tau;
    bool run_incoherent = sp.alpha > 0.000001;
    bool run_coherent = abs(sp.alpha - (PI / 4)) > 0.000001;

    for (int i = 0; i < sr.n_runs; i++) {
        sp.magnetic_field = range.x + step_size * i;
        sr.xs[i] = sp.magnetic_field;
        sr.xs_temperature[i] = sp.temperature;


        sp.mode == MODE_SIMULATION;

        {
            double total_sigma_xx = 0;
            double total_sigma_xxi = 0;
            double total_sigma_xx_sq = 0;

            for (int j = 0; j < sr.iterations_per_run; j++) {
                sp.impurity_seed = 1123 + j * 831;
                double result, resulti;

                if (run_incoherent) {
                    sp.is_incoherent = 1;
                    es.Compute(sp, resulti);
                }
                else {
                    resulti = 0;
                }

                if (run_coherent) {
                    sp.is_incoherent = 0;
                    sp.tau = coherent_tau;
                    es.Compute(sp, result);
                }
                else {
                    result = 0;
                }

                double sxx = result / 1e8;
                double sxx_i = resulti / 1e8;

                total_sigma_xx += sxx;
                total_sigma_xxi += sxx_i;
                total_sigma_xx_sq += (sxx_i + sxx) * (sxx_i + sxx);
            }

            {
                sr.results_xxi[i] = (total_sigma_xxi) / (double)(sr.iterations_per_run);
                sr.results_xx[i] = (total_sigma_xx) / (double)(sr.iterations_per_run);


                double sxx_sq_exp = total_sigma_xx_sq / (double)(sr.iterations_per_run) + 1e-15;
                double sxx_exp = (total_sigma_xx + total_sigma_xxi) / (double)(sr.iterations_per_run);

                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(sr.iterations_per_run - 1));
                sr.delta_xxi[i] = sxx_std / sxx_exp;
            }
        }

        // SIGMA XY
        {
            sp.mode = MODE_SIGMA_XY;

            double total_sigma_xy = 0;
            double total_sigma_xyi = 0;
            double total_sigma_xy_sq = 0;

            for (int j = 0; j < sr.iterations_per_run; j++) {
                sp.impurity_seed = 1123 + j * 831;
                double result, resulti;

                if (run_incoherent) {
                    sp.is_incoherent = 1;
                    es.Compute(sp, resulti);
                }
                else {
                    resulti = 0;
                }

                if (run_coherent) {
                    sp.is_incoherent = 0;
                    sp.tau = coherent_tau;
                    es.Compute(sp, result);
                }
                else {
                    result = 0;
                }

                double sxy = result / 1e8;
                double sxy_i = resulti / 1e8;

                total_sigma_xy += sxy;
                total_sigma_xyi += sxy_i;
            }

            {
                sr.results_xy[i] = total_sigma_xy / (double)sr.iterations_per_run;
                sr.results_xyi[i] = total_sigma_xyi / (double)sr.iterations_per_run;
            }
        }

        printf("Progress %d/%d\n", i, sr.n_runs);
    }
}

void ImGuiRender(ElasticScattering& es) {
    SimulationParameters old_sp = sp;

#if 1
    sp.impurity_density = 5.34e14;
    sp.impurity_radius = 1.11e-8;
    sp.region_extends = 1e-6;
    sp.region_size = 4e-6;
    sp.alpha = 0.3;
    sp.dim = 128;
    sp.tau = 1e-11;
#else
    sp.impurity_density = 5.34e12;
    sp.impurity_radius = 1.11e-7;
    sp.region_extends = 1e-5;
    sp.region_size = 4e-5;
    sp.alpha = PI / 4;
    sp.dim = 128;
#endif

    const std::vector<double> zs{ 15, 60 };

    for (int i = 0; i < zs.size(); i++) {
        sp.temperature = zs[i];

        QueryPerformanceCounter(&beginClock);

        SimulationResult result;
        ComputeSimulation(es, result);
        QueryPerformanceCounter(&endClock);
        result.time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        logger.LogResult(sp, result);
        printf("Simulation completed %d/%d\n", i + 1, zs.size());
    }

    sp = old_sp;
}

int app_main(const InitParameters& init)
{
    QueryPerformanceFrequency(&clockFrequency);

    sp = ParametersFactory::GenerateMinimal();
    GPUElasticScattering es(init);

    return 0;
}

