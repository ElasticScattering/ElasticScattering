#include "SimulationRunner.h"

#include "Logger.h"
#include "SimulationResult.h"
#include "scattering/Simulation.h"
#include "scattering/Grid.h"
#include "scattering/escl/constants.h"

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

void SimulationRunner::Run(const InitParameters& init)
{
    ParseConfig(init.config_file);
    CreateOutputDirectories();
    PrintSimulationInfo();

    QueryPerformanceFrequency(&clockFrequency);

    auto sp = cfg.scattering_params;
    ScatteringParameters sp_inc, sp_coh;
    
    sp_inc = cfg.scattering_params;
    sp_inc.is_coherent = false;
    CompleteSimulationParameters(sp_inc);

    sp_coh = cfg.scattering_params;
    sp_coh.is_coherent = true;
    CompleteSimulationParameters(sp_coh);
    
    const int S = cfg.num_samples;

    std::vector<SampleResult> sample_results_coh(S);
    std::vector<SampleResult> sample_results_inc(S);

    std::random_device random_device;
    SimulationCPU es;

    for (int i = 0; i < S; i++) {
        QueryPerformanceCounter(&beginClock);
        auto grid = Grid(sp.impurity_count, random_device(), sp.impurity_spawn_range, sp.impurity_radius, sp.cells_per_row);
        QueryPerformanceCounter(&endClock);

        GlobalMetrics gm;
        gm.coherent = sp.is_coherent;
        gm.particles_per_row                   = sp.dim - 1;
        gm.cells_per_row                       = sp.cells_per_row;
        gm.phi_values                          = 4 * sp.integrand_steps;
        gm.unique_impurity_count               = sp.impurity_count;
        gm.additional_impurities               = grid.GetImpurities().size() - gm.unique_impurity_count;
        gm.grid_time_elapsed                   = GetElapsedTime();
        gm.avg_impurities_in_cell              = grid.GetImpurities().size() / (double)(sp.cells_per_row * sp.cells_per_row);
        gm.avg_impuritiies_in_cell_overlapping = gm.avg_impurities_in_cell - sp.max_expected_impurities_in_cell;

        Logger::CreateMetricsLog(GetMetricsPath(i, gm.coherent), gm);

        printf(" ");

        auto sample_coh = RunSample(i, es, sp_coh, grid);
        sample_results_coh[i] = sample_coh;

        auto sample_inc = RunSample(i, es, sp_inc, grid);
        sample_results_inc[i] = sample_inc;

        printf("\n");
    }

    // Sommeer alle samples en log resultaat.
    for (int i = 0; i < cfg.temperatures.size(); i++) {
        for (int j = 0; j < cfg.magnetic_field_range.n; j++) {
            Sigma coherent, incoherent;
            double dxx_squared = 0;

            for (int s = 0; s < S; s++) {
                // Sommeer <MF,T> voor alle samples.
                auto coh = sample_results_coh[s].results[i][j];
                auto inc = sample_results_inc[s].results[i][j];

                coherent   += coh;
                incoherent += inc;

                dxx_squared += pow(coh.xx + inc.xx, 2);
            }

            coherent.xx   /= S;
            coherent.xy   /= S;
            incoherent.xx /= S;
            incoherent.xy /= S;

            DataRow row;
            row.temperature    = cfg.temperatures[i];
            row.magnetic_field = cfg.magnetic_field_range.min + j * cfg.magnetic_field_range.step_size;
            row.coherent       = coherent;
            row.coherent       = incoherent;

            double sxx_sq_exp = dxx_squared / S + 1e-15;
            double sxx_exp    = (coherent.xx + incoherent.xx) / S;
            double sxx_std    = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (S - 1));

            row.xxd           = sxx_std / sxx_exp;
            Logger::LogResult(GetResultPath(i), row);
        }
    }
}

SampleResult SimulationRunner::RunSample(const int sample_index, Simulation& es, ScatteringParameters& sp, const Grid& grid)
{
    const int T = cfg.temperatures.size();
    const int N = cfg.magnetic_field_range.n;
    SampleResult sr(T, N);

    auto metrics_path = GetMetricsPath(sample_index, sp.is_coherent);

    for (int j = 0; j < N; j++) {
        UpdateMagneticField(sp, cfg.magnetic_field_range.min + j * cfg.magnetic_field_range.step_size);
        Metrics metrics(j);

        QueryPerformanceCounter(&beginClock);
        es.ComputeLifetimes(sp, grid, metrics);
        QueryPerformanceCounter(&endClock);
        metrics.time_elapsed_lifetimes = GetElapsedTime();
        
        QueryPerformanceCounter(&beginClock);
        for (int i = 0; i < T; i++) {
            UpdateTemperature(sp, cfg.temperatures[i]);
            auto iteration = es.DeriveTemperature(sp.temperature);

            sr.results[i][j] = iteration.result;
            auto image_path = GetImagePath(i, j, sample_index, sp.is_coherent);
            Logger::LogImages(image_path, sp.dim - 1, iteration);
        }
        QueryPerformanceCounter(&endClock);
        metrics.time_elapsed_temperatures = GetElapsedTime();
        
        double n_lt = pow(sp.dim - 1, 2) * 4.0 * sp.integrand_steps;
        metrics.pct_particles_inside_impurity  = (double)(100 * metrics.particles_inside_impurity) / n_lt;
        metrics.prt_cells_passed               = (double)metrics.cells_passed / n_lt;
        metrics.pct_prt_cells_passed           = 100 * metrics.prt_cells_passed / pow(sp.cells_per_row, 2);
        metrics.prt_impurity_intersections     = (double)metrics.impurity_intersections / n_lt;
        metrics.pct_prt_impurity_intersections = 100 * metrics.prt_impurity_intersections / sp.impurity_count;
        metrics.mln_cells_passed               = (double)metrics.cells_passed / 1'000'000.0;
        metrics.mln_impurity_intersections     = (double)metrics.impurity_intersections / 1'000'000.0;

        Logger::LogMetrics(metrics_path, metrics);
        printf(".");
    }

    return sr;
}

void SimulationRunner::UpdateMagneticField(ScatteringParameters& sp, const double magnetic_field) {
    sp.magnetic_field = magnetic_field;
    sp.angular_speed = E * sp.magnetic_field / M;
}

void SimulationRunner::UpdateTemperature(ScatteringParameters& sp, const double temperature) {
    sp.temperature = temperature;

    if (!sp.is_coherent) {
        sp.tau = HBAR / (KB * sp.temperature);
        sp.default_max_lifetime = 15.0 * sp.tau;
    }
}

void SimulationRunner::CompleteSimulationParameters(ScatteringParameters& sp) {
    sp.angular_speed = E * sp.magnetic_field / M;

    if (!sp.is_coherent) sp.tau = HBAR / (KB * sp.temperature);
    sp.default_max_lifetime = 15.0 * sp.tau;

    {
        sp.impurity_spawn_range = { -sp.region_extends, sp.region_size + sp.region_extends };
        double area_length = sp.impurity_spawn_range.y - sp.impurity_spawn_range.x;
        sp.impurity_count = max(1, (int)ceil(area_length * area_length * sp.impurity_density));
        sp.cells_per_row = max((int)ceil(sqrt(sp.impurity_count / (double)sp.max_expected_impurities_in_cell)), 1);
        sp.cell_size = (sp.impurity_spawn_range.y - sp.impurity_spawn_range.x) / (double)sp.cells_per_row;
    }

    {
        const double incoherent_area = sp.alpha * 2.0;
        sp.integrand_angle_area = !sp.is_coherent ? incoherent_area : (PI / 2.0 - incoherent_area);
        sp.integrand_step_size = sp.integrand_angle_area / (sp.integrand_steps - 1);

        sp.integrand_start_angle = (!sp.is_coherent ? -sp.alpha : sp.alpha);
    }
}

void SimulationRunner::PrintSimulationInfo() const
{
    printf("Starting simulation (%i samples)\n", cfg.num_samples);
    printf("\tWriting to: %s\n\n", cfg.base_output_directory.c_str());


    printf("/%s\\\n", std::string(cfg.magnetic_field_range.n * 2, '-').c_str());
}

std::string SimulationRunner::GetAvailableDirectory(std::string base)
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

void SimulationRunner::ParseConfig(std::string file)
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

    Range magnetic_field;
    magnetic_field.min = atof(values.at("magnetic_field_min").c_str());
    magnetic_field.max = atof(values.at("magnetic_field_max").c_str());
    magnetic_field.n = atoi(values.at("magnetic_field_n").c_str());
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

    cfg.num_samples = atoi(values.at("num_samples").c_str());
    cfg.base_output_directory = GetAvailableDirectory(values.at("output_directory"));

    auto sp = cfg.scattering_params;
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

    sp.is_clockwise = values.at("clockwise").c_str() == "true";

    CompleteSimulationParameters(sp);

    cfg.scattering_params = sp;
}

void SimulationRunner::CreateOutputDirectories() const
{
    for (int i = 0; i < cfg.num_samples; i++)
    {
        auto sample_path = GetSamplePath(i);
        std::filesystem::create_directory(sample_path);
        std::filesystem::create_directory(sample_path + "/Incoherent/");
        std::filesystem::create_directory(sample_path + "/Coherent/");
    }

    for (int i = 0; i < cfg.temperatures.size(); i++)
        Logger::CreateResultLog(GetResultPath(i), cfg, cfg.temperatures[i]);
}

std::string SimulationRunner::GetSamplePath(int sample_idx) const
{
    return cfg.base_output_directory + "/Sample " + std::to_string(sample_idx);
}

std::string SimulationRunner::GetImagePath(int t_idx, int m_idx, int sample_idx, bool coherent) const
{
    auto type = coherent ? "/Coherent/" : "/Incoherent/";
    return GetSamplePath(sample_idx) + type + "T" + std::to_string(t_idx) + " MF" + std::to_string(m_idx) + ".png";
}

std::string SimulationRunner::GetMetricsPath(int sample_idx, bool coherent) const
{
    auto type = coherent ? "Coherent" : "Incoherent";
    return GetSamplePath(sample_idx) + "/Metrics " + type + ".txt";
}

std::string SimulationRunner::GetResultPath(int t_idx) const
{
    return cfg.base_output_directory + "/T" + std::to_string(t_idx) + ".dat";
}
