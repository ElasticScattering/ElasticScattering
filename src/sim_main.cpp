
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "sim_main.h"

#include "SimulationConfiguration.h"
#include "SimulationResult.h"
#include "utils/ParametersFactory.h"
#include "scattering/escl/constants.h"

//#include <stb_image_write/stb_image_write.h>
#include <random>
#include <windows.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <iomanip> 
#include <ctime>

#include <unordered_map>

int sim_main(const InitParameters& init)
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    
    SimulationConfiguration cfg = ParseConfig("default.config");

    const std::vector<double> temperatures { 15, 60 };

    ElasticScatteringCPU es;

    for (int i = 0; i < temperatures.size(); i++) {
        cfg.scattering_params.temperature = temperatures[i];

        QueryPerformanceCounter(&beginClock);
        auto simulation_result = RunSimulation(es, cfg);
        QueryPerformanceCounter(&endClock);
        simulation_result.time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        LogResult(cfg, simulation_result);
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

    SimulationResult sr(sp.number_of_runs);
    
    std::random_device random_device;
    const int base_seed = random_device();

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

SimulationConfiguration& ParseConfig(std::string file)
{
    std::unordered_map<std::string, std::string> values;

    std::filebuf fb;
    if (!fb.open(file, std::ios::in)) {
        printf("Could not load config file.");
        exit(-1);
    }

    std::istream is_file(&fb);

    std::string line;
    while (std::getline(is_file, line))
    {
        std::istringstream is_line(line);
        std::string line;
        if (std::getline(is_line, line))
        {
            if (line.size() == 0 || line[0] == '\n' || line[0] == '#' || (line.size() > 1 && line[0] == ':' && line[1] == '/'))
                continue;


            auto pos = line.find(' ', 0);
            auto pos2 = line.rfind(' ', line.size());

            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos2+1, line.size()-1);;

            values.insert(std::unordered_map<std::string, std::string>::value_type(key, value));
        }
    }

    fb.close();

    SimulationConfiguration cfg;

    cfg.magnetic_field_min = atof(values.at("magnetic_field_min").c_str());
    cfg.magnetic_field_max = atof(values.at("magnetic_field_max").c_str());
    cfg.number_of_runs = atoi(values.at("number_of_runs").c_str());
    cfg.samples_per_run = atoi(values.at("samples_per_run").c_str());

    ScatteringParameters sp;
    sp.integrand_steps = atoi(values.at("integrand_steps").c_str());
    sp.dim = atoi(values.at("dimension").c_str());;

    //sp.temperature = atof(values.at("temperature").c_str());
    sp.tau = atof(values.at("tau").c_str());
    //sp.magnetic_field = atof(values.at("magnetic_field").c_str());
    sp.alpha = atof(values.at("alpha").c_str());
    sp.particle_speed = atof(values.at("particle_speed").c_str());

    sp.impurity_density = atof(values.at("impurity_density").c_str());
    sp.impurity_radius = atof(values.at("impurity_radius").c_str());
    sp.region_extends = atof(values.at("region_extends").c_str());
    sp.region_size = atof(values.at("region_size").c_str());
    sp.max_expected_impurities_in_cell = atoi(values.at("max_expected_impurities_in_cell").c_str());

    sp.is_diag_regions = atoi(values.at("is_diag_regions").c_str());
    sp.is_clockwise = atoi(values.at("is_clockwise").c_str());
    //sp.is_incoherent = atoi(values.at("is_incoherent").c_str());

    cfg.scattering_params = sp;

    return cfg;
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
