#include "sim_main.h"

#include "scattering/Simulation.h"
#include "scattering/escl/constants.h"
#include "Logger.h"

#include <string>
#include <unordered_map>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <iomanip> 

#include <random>
#include <windows.h>

void sim_main(const InitParameters& init)
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);
    
    SimulationConfiguration cfg = ParseConfig(init.config_file);
    PrintSimulationInfo(cfg);

    auto base_name = cfg.output_directory + "/results_";
    for (int i = 0; i < cfg.temperatures.size(); i++) 
        Logger::CreateLog(base_name + std::to_string(i) + ".dat", cfg, cfg.temperatures[i]);

    SimulationCPU es;
    for (int i = 0; i < cfg.magnetic_field_range.n; i++) {
        cfg.scattering_params.magnetic_field = cfg.magnetic_field_range.min + cfg.magnetic_field_range.step_size * i;

        cfg.intermediates_directory = cfg.output_directory + "/" + std::to_string(i) + ". Intermediates";
        std::filesystem::create_directory(cfg.intermediates_directory);

        QueryPerformanceCounter(&beginClock);
        RunSimulation(cfg, es);
        QueryPerformanceCounter(&endClock);
        double time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        printf("Simulation completed (%d/%d)\n", i + 1, cfg.magnetic_field_range.n);
    }
}

void RunSimulation(const SimulationConfiguration& cfg, Simulation& es)
{
    ScatteringParameters sp_inc = cfg.scattering_params;
    sp_inc.is_incoherent = 1;

    ScatteringParameters sp_coh = cfg.scattering_params;
    sp_coh.is_incoherent = 0;

    CompleteSimulationParameters(sp_inc);
    CompleteSimulationParameters(sp_coh);

    std::random_device random_device;

    std::vector<Sigma>  coherent_results(cfg.temperatures.size());
    std::vector<Sigma>  incoherent_results(cfg.temperatures.size());
    std::vector<double> total_dxx_sq(cfg.temperatures.size());

    for (int j = 0; j < cfg.samples_per_run; j++) {
        auto impurity_index = Grid(sp_coh.impurity_count, random_device(), sp_coh.impurity_spawn_range, sp_coh.impurity_radius, sp_coh.cells_per_row);

        std::vector<Sigma> sample_coherent_results(cfg.temperatures.size());
        std::vector<Sigma> sample_incoherent_results(cfg.temperatures.size());

        es.ComputeLifetimes(sp_coh, impurity_index);
        for (int i = 0; i < cfg.temperatures.size(); i++) {
            UpdateTemperature(sp_coh, cfg.temperatures[i]);
            auto iteration = es.DeriveTemperature(cfg.temperatures[i]);

            coherent_results[i]       += iteration.result;
            sample_coherent_results[i] = iteration.result;

            Logger::LogImages(cfg.intermediates_directory + "/" + std::to_string(i) + ". Sample " + std::to_string(j) + ".png", sp_coh.dim - 1, sp_inc.tau, iteration);
        }

        es.ComputeLifetimes(sp_inc, impurity_index);
        for (int i = 0; i < cfg.temperatures.size(); i++) {
            UpdateTemperature(sp_inc, cfg.temperatures[i]);
            auto iteration = es.DeriveTemperature(cfg.temperatures[i]);

            incoherent_results[i]       += iteration.result;
            sample_incoherent_results[i] = iteration.result;
            Logger::LogImages(cfg.intermediates_directory + "/" + std::to_string(i) + ". Sample " + std::to_string(j) + "I.png", sp_inc.dim - 1, sp_inc.tau, iteration);
        }

        for (int i = 1; i < cfg.temperatures.size(); i++) {
            total_dxx_sq[i] += pow(sample_coherent_results[i].xx + sample_incoherent_results[i].xx, 2);
        }

        printf("Sample completed!\n");
    }

    double sample_count = (double)cfg.samples_per_run;
    for (int i = 0; i < cfg.temperatures.size(); i++)
    {
        coherent_results[i].xx /= sample_count;
        coherent_results[i].xy /= sample_count;

        incoherent_results[i].xx /= sample_count;
        incoherent_results[i].xy /= sample_count;
    }

    for (int i = 0; i < cfg.temperatures.size(); i++)
    {
        DataRow row;
        row.magnetic_field = cfg.scattering_params.magnetic_field;
        row.temperature    = cfg.temperatures[i];
        row.coherent       = coherent_results[i];
        row.incoherent     = incoherent_results[i];

        double sxx_sq_exp = total_dxx_sq[i] / sample_count + 1e-15;
        double sxx_exp    = (coherent_results[i].xx + incoherent_results[i].xx) / sample_count;
        double sxx_std    = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (sample_count - 1));

        row.xxd = sxx_std / sxx_exp;
        
        Logger::LogResult(cfg.output_directory + "/results_" + std::to_string(i) + ".dat", row);
    }
}

void UpdateTemperature(ScatteringParameters& sp, double temperature) {
    sp.temperature = temperature;

    if (sp.is_incoherent == 1) {
        sp.tau = HBAR / (KB * sp.temperature);
        sp.default_max_lifetime = 15.0 * sp.tau;
    }
}

void CompleteSimulationParameters(ScatteringParameters& sp) {
    sp.angular_speed = E * sp.magnetic_field / M;

    if (sp.is_incoherent == 1) sp.tau = HBAR / (KB * sp.temperature);
    sp.default_max_lifetime = 15.0 * sp.tau;

    {
        sp.impurity_spawn_range = { -sp.region_extends, sp.region_size + sp.region_extends };
        double area_length = sp.impurity_spawn_range.y - sp.impurity_spawn_range.x;
        sp.impurity_count = max(1, (int)ceil(area_length * area_length * sp.impurity_density));
        sp.cells_per_row = max((int)ceil(sqrt(sp.impurity_count / (double)sp.max_expected_impurities_in_cell)), 1);
        sp.cell_size = (sp.impurity_spawn_range.y - sp.impurity_spawn_range.x) / (double)sp.cells_per_row;
    }

    {
        bool incoherent = (sp.is_incoherent == 1);

        const double incoherent_area = sp.alpha * 2.0;
        sp.integrand_angle_area = incoherent ? incoherent_area : (PI / 2.0 - incoherent_area);
        sp.integrand_step_size = sp.integrand_angle_area / (sp.integrand_steps - 1);

        sp.integrand_start_angle = (incoherent ? -sp.alpha : sp.alpha);
    }
}

std::string GetAvailableDirectory(std::string base)
{
    std::filesystem::create_directory(base);
    std::string base_dir_name = base + "/Result_";

    unsigned int n = 0;
    while (true)
    {
        std::string dir = base_dir_name + std::to_string(n);
        if (!std::filesystem::exists(dir))
        {
            std::filesystem::create_directory(dir);
            return dir;
        }

        n++;
    }
}

SimulationConfiguration ParseConfig(std::string file)
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
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1, line.size() - 1);

            values.insert(std::unordered_map<std::string, std::string>::value_type(key, value));
        }
    }

    fb.close();

    SimulationConfiguration cfg;

    Range magnetic_field;
    magnetic_field.min       = atof(values.at("magnetic_field_min").c_str());
    magnetic_field.max       = atof(values.at("magnetic_field_max").c_str());
    magnetic_field.n         = atoi(values.at("magnetic_field_n").c_str());
    magnetic_field.step_size = (magnetic_field.max - magnetic_field.min) / (double)(magnetic_field.n - 1);
    cfg.magnetic_field_range = magnetic_field;

    {
        std::stringstream string_stream(values.at("temperatures"));

        int i = 0;
        while (string_stream.good())
        {
            std::string a;
            std::getline(string_stream, a, ' ');
            cfg.temperatures.push_back(atof(a.c_str()));
            i++;
        }
    }

    cfg.samples_per_run  = atoi(values.at("samples_per_run").c_str());
    cfg.output_directory = GetAvailableDirectory(values.at("output_directory"));

    ScatteringParameters sp;
    sp.integrand_steps                 = atoi(values.at("integrand_steps").c_str());
    sp.dim                             = atoi(values.at("dimension").c_str());;
    sp.values_per_particle             = 4 * sp.integrand_steps;

    sp.tau                             = atof(values.at("tau").c_str());
    sp.alpha                           = atof(values.at("alpha").c_str());
    sp.particle_speed                  = atof(values.at("particle_speed").c_str());

    sp.impurity_density                = atof(values.at("impurity_density").c_str());
    sp.impurity_radius                 = atof(values.at("impurity_radius").c_str());
    sp.region_extends                  = atof(values.at("region_extends").c_str());
    sp.region_size                     = atof(values.at("region_size").c_str());
    sp.max_expected_impurities_in_cell = atoi(values.at("max_expected_impurities_in_cell").c_str());

    sp.is_clockwise                    = atoi(values.at("is_clockwise").c_str());

    cfg.scattering_params = sp;

    return cfg;
}

void PrintSimulationInfo(const SimulationConfiguration &cfg)
{
    printf("Starting simulation...\n");
    printf("\tWriting to: %s\n\n", cfg.output_directory.c_str());
}
