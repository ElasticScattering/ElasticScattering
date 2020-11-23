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

    SimulationConfiguration cfg = ParseConfig("default.config"); // Naam via cmd?
    printf("Writing to: %s\n", cfg.output_directory.c_str());

    const std::vector<double> temperatures { 15, 60 }; // @Todo, dit ook inlezen?

    for (int i = 0; i < temperatures.size(); i++) {
        cfg.scattering_params.temperature = temperatures[i];
        
        CreateLog(cfg);

        QueryPerformanceCounter(&beginClock);
        RunSimulation(cfg);
        QueryPerformanceCounter(&endClock);
        
        double time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
        FinishLog(cfg.result_file, time_elapsed);
       
        printf("Simulation completed (%d/%d)\n", i + 1, temperatures.size());
    }
}

IterationResult RunIteration(const ScatteringParameters& sp, const ImpurityIndex& grid)
{
    auto lifetimes          = ElasticScatteringCPU::ComputeLifetimes(sp, grid);
    auto particle_lifetimes = ElasticScatteringCPU::IntegrateParticle(sp, lifetimes);
    auto particle_sigmas    = ElasticScatteringCPU::ComputeSigmas(sp, lifetimes);

    auto result             = ElasticScatteringCPU::IntegrateResult(sp, lifetimes);
    result.xx /= 1e8;
    result.xy /= 1e8;

    IterationResult b;
    b.result = result;
    b.lifetimes = particle_lifetimes;
    b.sigma = particle_sigmas;
    
    return b;
}

void RunSimulation(SimulationConfiguration& cfg)
{
    ScatteringParameters sp = cfg.scattering_params;
    double coherent_tau = sp.tau;
    bool run_incoherent = sp.alpha > 0.000001;
    bool run_coherent = abs(sp.alpha - (PI / 4)) > 0.000001;

    std::random_device random_device;

    double mf_step_size = (cfg.magnetic_field_max - cfg.magnetic_field_min) / cfg.number_of_runs;
    for (int i = 0; i < cfg.number_of_runs; i++) {
        ElasticScattering::UpdateSimulationParameters(sp, cfg.magnetic_field_min + mf_step_size * i, sp.temperature);

        {
            SigmaResult total_coherent, total_incoherent;
            double total_sigma_xx_sq = 0;

            for (int j = 0; j < cfg.samples_per_run; j++) {
                auto impurity_index = ImpurityIndex(sp.impurity_count, random_device(), sp.impurity_spawn_range, sp.impurity_radius, sp.cells_per_row);

                sp.is_incoherent = 1;
                auto incoherent = RunIteration(sp, impurity_index);
                
                sp.is_incoherent = 0;
                sp.tau = coherent_tau;
                auto coherent = RunIteration(sp, impurity_index);

                LogImages(cfg.output_directory + "/" + std::to_string(i) + ". Sample " + std::to_string(j) + ".png", sp.dim-1, sp.tau, coherent, incoherent);

                /*
                if (run_incoherent) {
                    sp.is_incoherent = 1;
                    incoherent = RunIteration(sp, impurity_index);
                }

                if (run_coherent) {
                    sp.is_incoherent = 0;
                    sp.tau = coherent_tau;
                    coherent = RunIteration(sp, impurity_index);
                }
                */

                total_coherent.xx += coherent.result.xx;
                total_coherent.xy += coherent.result.xy;
                total_incoherent.xx += incoherent.result.xx;
                total_incoherent.xy += incoherent.result.xy;

                total_sigma_xx_sq += (incoherent.result.xx + coherent.result.xx) * (incoherent.result.xx + coherent.result.xx);
            }

            DataRow row(total_coherent, total_incoherent, cfg.samples_per_run);
            row.temperature = sp.temperature;
            row.magnetic_field = sp.magnetic_field;

            {
                double sxx_sq_exp = total_sigma_xx_sq / (double)(cfg.samples_per_run) + 1e-15;
                double sxx_exp = (total_coherent.xx + total_incoherent.xx) / (double)(cfg.samples_per_run);

                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(cfg.samples_per_run - 1));
                
                row.xxd = sxx_std / sxx_exp;
            }
         
            LogResult(cfg.result_file, row);
        }
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
    cfg.output_directory = GetAvailableDirectory(values.at("output_directory"));
    cfg.result_file = cfg.output_directory + std::string("/_results.dat");

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
    
    ElasticScattering::CompleteSimulationParameters(sp);

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

void CreateLog(const SimulationConfiguration& cfg)
{
    std::ofstream file;
    file.open(cfg.result_file);

    file << "# Elastic Scattering simulation results" << std::endl;
    file << "# Number of runs: " << cfg.number_of_runs << std::endl;
    file << "# Samples per row: " << cfg.samples_per_run << std::endl;
    
    ScatteringParameters sp = cfg.scattering_params;

    file << std::scientific << std::setprecision(3);
    file << "#" << std::endl;
    file << "###########################" << std::endl;
    file << "# Scattering parameters:  #" << std::endl;
    file << "###########################" << std::endl;
    file << "#\t" << "Dimension:        " << sp.dim << std::endl;
    file << "#\t" << "Integrand steps:  " << sp.integrand_steps << std::endl;
    file << "#\t" << "Diag. regions:    " << ((sp.is_diag_regions == 1) ? "True" : "False") << std::endl;
    file << "#\t" << "Clockwise:        " << ((sp.is_clockwise == 1) ? "True" : "False") << std::endl;
    file << "#" << std::endl;
    //file << "#\t" << "Temperature:      " << sp.temperature << std::endl;
    file << "#\t" << "Tau:              " << sp.tau << std::endl;
    //file << "#\t" << "Magnetic field:   " << sp.magnetic_field << std::endl;
    file << "#\t" << "Alpha:            " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed:   " << sp.particle_speed << std::endl;
    //file << "#\t" << "Angular speed:    " << sp.angular_speed << std::endl;
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

void LogResult(const std::string file_path, const DataRow& row)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);
    
    file << std::scientific << std::setprecision(10);
    file << row.temperature << " " << row.incoherent.xx << " " << row.coherent.xx << " " << row.incoherent.xy << " " << row.coherent.xy << " " << row.xxd << std::endl;
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

void LogImages(const std::string file, const int dim, const double scale, const IterationResult coherent, const IterationResult incoherent)
{
    const int image_width  = dim * 3;
    const int image_height = dim * 2;
    const int dim2 = dim * 2;
    
    std::vector<unsigned char> pixels(dim * image_width * image_height);

    int i, j = 0;

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + i] = (unsigned char)(coherent.lifetimes[j * dim + i] / scale * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + (i + dim)] = (unsigned char) (coherent.sigma.sigma_xx[j * dim + i] / 3.0 * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + (i + dim2)] = (unsigned char)(coherent.sigma.sigma_xy[j * dim + i] / 1.5 * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[(j + dim) * image_width + i] = (unsigned char)(205 - incoherent.lifetimes[j * dim + i] / scale);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[(j + dim) * image_width + (i + dim)] = (unsigned char)(incoherent.sigma.sigma_xx[j * dim + i] * 3.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[(j + dim) * image_width + (i + dim2)] = (unsigned char)(incoherent.sigma.sigma_xy[j * dim + i] * 1.5);

    stbi_write_png(file.c_str(), image_width, image_height, 1, pixels.data(), image_width);
}

/*
void LogImage(const std::string file, const int dim, const double scale, const std::vector<double> data)
{
    std::vector<unsigned char> pixels(dim * dim * 3);
    int pix_idx = 0;
    for (int j = 0; j < dim; j++)
        for (int i = 0; i < dim; i++)
        {
            auto k = (unsigned char)(data[j * dim + i] / scale);
            pixels[pix_idx++] = k;
        }

    stbi_write_png(file.c_str(), dim, dim, 3, pixels.data(), 1);
}*/