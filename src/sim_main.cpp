#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "sim_main.h"

#include "scattering/ElasticScattering.h"
#include "scattering/escl/constants.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

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

void sim_main(const InitParameters& init)
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);
    
    SimulationConfiguration cfg = ParseConfig(init.config_file);
    printf("Writing to: %s\n", cfg.output_directory.c_str());

    auto base_name = cfg.output_directory + "/results_";
    for (int i = 0; i < cfg.temperatures.size(); i++) {
        cfg.result_file = base_name + std::to_string(i) + ".dat";
        CreateLog(cfg, cfg.temperatures[i]);
    }

    ElasticScatteringCPU es;
    for (int i = 0; i < cfg.magnetic_field_range.n; i++) {
        cfg.scattering_params.magnetic_field = cfg.magnetic_field_range.min + cfg.magnetic_field_range.step_size * i;

        cfg.intermediates_directory = cfg.output_directory + "/" + std::to_string(i) + ". Intermediates";
        std::filesystem::create_directory(cfg.intermediates_directory);

        QueryPerformanceCounter(&beginClock);
        RunSimulation(cfg, es);
        QueryPerformanceCounter(&endClock);
        double time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        //FinishLog(cfg.result_file, time_elapsed);
        printf("Simulation completed (%d/%d)\n", i + 1, cfg.magnetic_field_range.n);
    }
}

void RunSimulation(const SimulationConfiguration& cfg, ElasticScattering& es)
{

    ScatteringParameters sp_inc = cfg.scattering_params;
    sp_inc.is_incoherent = 1;

    ScatteringParameters sp_coh = cfg.scattering_params;
    sp_coh.is_incoherent = 0;

    es.CompleteSimulationParameters(sp_inc);
    es.CompleteSimulationParameters(sp_coh);

    std::random_device random_device;

    std::vector<Sigma>  coherent_results(cfg.temperatures.size());
    std::vector<Sigma>  incoherent_results(cfg.temperatures.size());
    std::vector<double> total_dxx_sq(cfg.temperatures.size());

    for (int j = 0; j < cfg.samples_per_run; j++) {
        auto impurity_index = ImpurityIndex(sp_coh.impurity_count, random_device(), sp_coh.impurity_spawn_range, sp_coh.impurity_radius, sp_coh.cells_per_row);

        std::vector<Sigma> sample_coherent_results(cfg.temperatures.size());
        std::vector<Sigma> sample_incoherent_results(cfg.temperatures.size());

        es.ComputeLifetimes(sp_coh, impurity_index);
        for (int i = 0; i < cfg.temperatures.size(); i++) {
            auto iteration = es.DeriveTemperature(cfg.temperatures[i]);

            coherent_results[i]       += iteration.result;
            sample_coherent_results[i] = iteration.result;

            LogImages(cfg.intermediates_directory + "/" + std::to_string(i) + ". Sample " + std::to_string(j) + ".png", sp_coh.dim - 1, sp_inc.tau, iteration);
        }

        es.ComputeLifetimes(sp_inc, impurity_index);
        for (int i = 0; i < cfg.temperatures.size(); i++) {
            auto iteration = es.DeriveTemperature(cfg.temperatures[i]);

            incoherent_results[i]       += iteration.result;
            sample_incoherent_results[i] = iteration.result;
            LogImages(cfg.intermediates_directory + "/" + std::to_string(i) + ". Sample " + std::to_string(j) + "I.png", sp_inc.dim - 1, sp_inc.tau, iteration);
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
        
        LogResult(cfg.output_directory + "/results_" + std::to_string(i) + ".dat", row);
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
            std::string value = line.substr(pos+1, line.size()-1);

            values.insert(std::unordered_map<std::string, std::string>::value_type(key, value));
        }
    }

    fb.close();

    SimulationConfiguration cfg;

    Range magnetic_field;
    magnetic_field.min = atof(values.at("magnetic_field_min").c_str());
    magnetic_field.max = atof(values.at("magnetic_field_max").c_str());
    magnetic_field.n = atoi(values.at("magnetic_field_n").c_str());
    magnetic_field.step_size = (magnetic_field.max - magnetic_field.min) / (double)(magnetic_field.n - 1);
    cfg.magnetic_field_range = magnetic_field;

    cfg.temperatures.push_back(atof(values.at("temperatures").c_str()));

    cfg.samples_per_run = atoi(values.at("samples_per_run").c_str());
    cfg.output_directory = GetAvailableDirectory(values.at("output_directory"));

    ScatteringParameters sp;
    sp.integrand_steps = atoi(values.at("integrand_steps").c_str());
    sp.dim = atoi(values.at("dimension").c_str());;
    sp.values_per_particle = 4 * sp.integrand_steps;

    sp.tau = atof(values.at("tau").c_str());
    sp.alpha = atof(values.at("alpha").c_str());
    sp.particle_speed = atof(values.at("particle_speed").c_str());

    sp.impurity_density = atof(values.at("impurity_density").c_str());
    sp.impurity_radius = atof(values.at("impurity_radius").c_str());
    sp.region_extends = atof(values.at("region_extends").c_str());
    sp.region_size = atof(values.at("region_size").c_str());
    sp.max_expected_impurities_in_cell = atoi(values.at("max_expected_impurities_in_cell").c_str());

    sp.is_clockwise = atoi(values.at("is_clockwise").c_str());
    
    cfg.scattering_params = sp;

    return cfg;
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

void CreateLog(const SimulationConfiguration& cfg, double temperature)
{
    std::ofstream file;
    file.open(cfg.result_file);

    file << "# Elastic Scattering simulation results" << std::endl;
    file << "# Number of runs: " << cfg.magnetic_field_range.n << std::endl;
    file << "# Samples per row: " << cfg.samples_per_run << std::endl;
    
    ScatteringParameters sp = cfg.scattering_params;

    file << std::scientific << std::setprecision(3);
    file << "#" << std::endl;
    file << "###########################" << std::endl;
    file << "# Scattering parameters:  #" << std::endl;
    file << "###########################" << std::endl;
    file << "#\t" << "Dimension:        " << sp.dim << std::endl;
    file << "#\t" << "Integrand steps:  " << sp.integrand_steps << std::endl;
    file << "#" << std::endl;
    file << "#\t" << "Temperature:      " << temperature << std::endl;
    file << "#\t" << "Tau:              " << sp.tau << std::endl;
    file << "#\t" << "Alpha:            " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed:   " << sp.particle_speed << std::endl;
    file << "#\t" << "Clockwise:        " << ((sp.is_clockwise == 1) ? "True" : "False") << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Count:            " << sp.impurity_count << std::endl;
    file << "#\t" << "Region size:      " << sp.region_size << std::endl;
    file << "#\t" << "Region extends:   " << sp.region_extends << std::endl;
    file << "#\t" << "Density:          " << sp.impurity_density << std::endl;
    file << "#\t" << "Radius:           " << sp.impurity_radius << std::endl;

    file << "#\n# Constants:" << std::endl;
    file << "#\t" << "Particle mass:    " << M << std::endl;
    file << "#\t" << "E:                " << E << std::endl;
    file << "#\t" << "HBAR:             " << HBAR << std::endl;
    file << "#\t" << "C:                " << C1 << std::endl;
    file << "#\t" << "KB:               " << KB << std::endl;
    file << "#" << std::endl;
    file << "###########################" << std::endl;
    file << "# Results:                #" << std::endl;
    file << "###########################" << std::endl;
    file << "magnetic_field sigma_xx_inc sigma_xx_coh sigma_xy_inc sigma_xy_coh delta_xx" << std::endl;
}

void LogResult(const std::string file_path, const DataRow row)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    file << std::scientific << std::setprecision(10);
    file << row.magnetic_field << " " << row.incoherent.xx << " " << row.coherent.xx << " " << row.incoherent.xy << " " << row.coherent.xy << " " << row.xxd << std::endl;
}

void FinishLog(const std::string file_path, const double time_elapsed)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    auto date_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    file << "\n\n\n###################################" << std::endl;
    file << "# Completed on: " << std::put_time(std::localtime(&date_time), "%F %T") << "." << std::endl;
    file << "# Elapsed time: " << time_elapsed << " seconds." << std::endl;
}

void LogImages(const std::string file, const int dim, const double scale, const IterationResult iteration)
{
    const int image_width  = dim * 3;
    const int image_height = dim;
    
    std::vector<unsigned char> pixels(image_width * image_height);

    int i, j = 0;
    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + i] = (unsigned char)(iteration.particle_lifetimes[j * dim + i] / scale * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + (i + dim)] = (unsigned char)(iteration.sigmas.xx_buffer[j * dim + i] / 3.0 * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + (i + dim * 2)] = (unsigned char)(iteration.sigmas.xy_buffer[j * dim + i] / 1.5 * 255.0);

    stbi_write_png(file.c_str(), image_width, image_height, 1, pixels.data(), image_width);
}
