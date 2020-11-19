#pragma once
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "sim_main.h"

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "utils/ParametersFactory.h"
#include "scattering/escl/constants.h"

#include <stb_image_write/stb_image_write.h>
#include <random>
#include <windows.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <iomanip> 
#include <ctime>

int sim_main(const InitParameters& init)
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    ScatteringParameters sp = ParametersFactory::GenerateSimulation();
    SimulationConfiguration sim_params;
    sim_params.number_of_runs = 10;
    sim_params.samples_per_run = 2;
    sim_params.magnetic_field_min = 0.01;
    sim_params.magnetic_field_max = 40;
    sim_params.scattering_params = sp;

    const std::vector<double> temperatures { 15, 60 };
    PrintInfo(sim_params, temperatures.size());

    ElasticScatteringCPU es;

    for (int i = 0; i < temperatures.size(); i++) {
        sp.temperature = temperatures[i];

        QueryPerformanceCounter(&beginClock);
        auto simulation_result = RunSimulation(es, sim_params);
        QueryPerformanceCounter(&endClock);
        simulation_result.time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        LogResult(sim_params, simulation_result);
        printf("Simulation completed (%d/%d)\n", i + 1, temperatures.size());
    }

    return 0;
}

SimulationResult& RunSimulation(ElasticScattering &es, SimulationConfiguration& sp)
{
    double coherent_tau = sp.scattering_params.tau; //@Refactor
    bool run_incoherent = sp.scattering_params.alpha > 0.000001;
    bool run_coherent = abs(sp.scattering_params.alpha - (PI / 4)) > 0.000001;
    
    double step_size = (sp.magnetic_field_max - sp.magnetic_field_min) / sp.number_of_runs;

    //printf("magnetic_field sigma_xx_inc sigma_xx_coh sigma_xy_inc sigma_xy_coh delta_xx\n");

    std::random_device random_device;

    SimulationResult sr(sp.number_of_runs);

    for (int i = 0; i < sp.number_of_runs; i++) {
        sp.scattering_params.magnetic_field = sp.magnetic_field_min + step_size * i;

        {
            double total_sigma_xx = 0;
            double total_sigma_xxi = 0;
            double total_sigma_xx_sq = 0;

            double total_sigma_xy = 0;
            double total_sigma_xyi = 0;
            double total_sigma_xy_sq = 0;

            for (int j = 0; j < sp.samples_per_run; j++) {
                sp.scattering_params.impurity_seed = random_device();
                SigmaResult result_coherent, result_incoherent;

                if (run_incoherent) {
                    sp.scattering_params.is_incoherent = 1;
                    result_incoherent = es.ComputeResult(sp.scattering_params);
                }

                if (run_coherent) {
                    sp.scattering_params.is_incoherent = 0;
                    sp.scattering_params.tau = coherent_tau;
                    result_coherent = es.ComputeResult(sp.scattering_params);
                }

                double sxx = result_coherent.xx / 1e8;
                double sxy = result_coherent.xy / 1e8;
                
                double sxx_i = result_incoherent.xx / 1e8;
                double sxy_i = result_incoherent.xy / 1e8;

                total_sigma_xx += sxx;
                total_sigma_xxi += sxx_i;
                total_sigma_xx_sq += (sxx_i + sxx) * (sxx_i + sxx);

                total_sigma_xy += sxy;
                total_sigma_xyi += sxy_i;
            }

            DataRow row;
            row.temperature = sp.scattering_params.temperature;
            row.magnetic_field = sp.scattering_params.magnetic_field;

            {
                double sxx_sq_exp = total_sigma_xx_sq / (double)(sp.samples_per_run) + 1e-15;
                double sxx_exp = (total_sigma_xx + total_sigma_xxi) / (double)(sp.samples_per_run);

                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(sp.samples_per_run - 1));
                
                row.coherent.xx = total_sigma_xx / (double)(sp.samples_per_run);
                row.incoherent.xx = total_sigma_xxi / (double)(sp.samples_per_run);

                row.xxd = sxx_std / sxx_exp;

                row.coherent.xy  = total_sigma_xy / (double)sp.samples_per_run;
                row.incoherent.xy = total_sigma_xyi / (double)sp.samples_per_run;
            }

            sr.results[i] = row;
            //printf("%f %e %e %e %e %f\n", row.temperature, row.incoherent.xx, row.coherent.xx, row.incoherent.xy, row.coherent.xy, row.xxd);
        }
    }

    return sr;
}

void PrintInfo(const SimulationConfiguration& sp, int count)
{
    const int imp_count = sp.scattering_params.impurity_density * pow(sp.scattering_params.region_extends + sp.scattering_params.region_size, 2);
    const long long intersects_each = sp.samples_per_run * imp_count * 2.0 * pow(sp.scattering_params.dim, 2) * sp.scattering_params.integrand_steps * 4.0;
    const long long intersects = sp.number_of_runs * intersects_each;
    printf("Simulation info:\n");
    printf("Total intersections: %e\n", intersects * (long long)count);
    printf("Time estimate per data point: %f minutes.\n", (float)intersects_each / 2e9 / 60);
    printf("Time estimate per temperature: %f hours \n", (float)intersects / 2e9 / 3600.0);
    printf("Time estimate total: %f hours \n\n", (float)intersects / 2e9 / 3600.0 * count);
}

void LogResult(const SimulationConfiguration& sim_params, const SimulationResult& sr)
{
    auto date_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::filesystem::create_directory("ESLogs");
    std::filesystem::create_directory("ESLogs");

    std::string base_file_name = "ESLogs/Result_";
    std::ofstream file;

    for (unsigned int n = 0; ; ++n) 
    {
        std::string dir = base_file_name + std::to_string(n);
        if (!std::filesystem::exists(dir))
        {
            std::filesystem::create_directory(dir);

            std::string fname = dir + std::string("results.dat");

            file.open(fname.c_str());
            break;
        }
    }

    ScatteringParameters sp = sim_params.scattering_params;

    file << std::scientific << std::setprecision(10);
    file << "# Elastic Scattering Results" << std::endl;
    file << "# Completed on: " << std::put_time(std::localtime(&date_time), "%F %T") << "." << std::endl;
    file << "# Elapsed time: " << sr.time_elapsed << " seconds." << std::endl;
    file << "# Each row is the average of " << sim_params.samples_per_run << " iterations with different impurity seeds): " << std::endl;
    file << "# Scattering parameters:" << std::endl;
    file << "#\t" << "Integrand steps:  " << sp.integrand_steps << std::endl;
    file << "#\t" << "Dimension:        " << sp.dim << std::endl;
    file << "#\t" << "Diag. regions:    " << ((sp.is_diag_regions == 1) ? "True" : "False") << std::endl;
    file << "#\t" << "Clockwise:        " << ((sp.is_clockwise == 1) ? "True" : "False") << std::endl;
    file << "#" << std::endl;
    file << "#\t" << "Temperature:      " << sp.temperature << std::endl;
    file << "#\t" << "Tau:              " << sp.tau << std::endl;
    file << "#\t" << "Magnetic field:   " << sp.magnetic_field << std::endl;
    file << "#\t" << "Alpha:            " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed:   " << sp.particle_speed << std::endl;
    file << "#\t" << "Angular speed:    " << sp.angular_speed << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Region size:      " << sp.region_size << std::endl;
    file << "#\t" << "Region extends:   " << sp.region_extends << std::endl;
    file << "#\t" << "Density:          " << sp.impurity_density << std::endl;
    file << "#\t" << "Count:            " << sp.impurity_count << std::endl;
    file << "#\t" << "Radius:           " << sp.impurity_radius << std::endl;

    file << "#\n#Constants:" << std::endl;
    file << "#\t" << "Particle mass: " << M << std::endl;
    file << "#\t" << "E:             " << E << std::endl;
    file << "#\t" << "HBAR:          " << HBAR << std::endl;
    file << "#\t" << "C:             " << C1 << std::endl;
    file << "#\t" << "KB:            " << KB << std::endl;

    file << "#\n#Results:\n" << std::endl;

    file << "magnetic_field sigma_xx_inc sigma_xx_coh sigma_xy_inc sigma_xy_coh delta_xx" << std::endl;

    int n = sr.results.size();

    int idx = 0;
    for (int i = 0; i < n; i++) {
        const auto row = sr.results[i];
        file << row.temperature << " " << row.incoherent.xx << " " << row.coherent.xx << " " << row.incoherent.xy << " " << row.coherent.xy << " " << row.xxd << std::endl;
    }

    file.flush();
    file.close();
}
