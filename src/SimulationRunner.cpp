#include "SimulationRunner.h"

#include "scattering/Simulation.h"
#include "scattering/Grid.h"
#include "SimulationResult.h"
#include "Logger.h"

#include <string>
#include <iostream>

#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>

#include <random>

void SimulationRunner::Run(const InitParameters& init)
{
    QueryPerformanceFrequency(&clockFrequency);

    ParseConfig(init.config_file);
    CreateOutputDirectories();
    PrintSimulationInfo();

    const int S = cfg.num_samples;

    std::vector<SampleResult> sample_results_coh(S);
    std::vector<SampleResult> sample_results_inc(S);

    std::random_device random_device;
    SimulationCPU es(cfg.particles_per_row-1, cfg.quadrant_integral_steps);

    const Settings& ss = cfg.settings;

    for (int i = 0; i < S; i++) {
        std::cout << " ";
        
        QueryPerformanceCounter(&beginClock);
        auto grid = Grid(random_device(), ss.region_size, ss.region_extends, ss.impurity_density, ss.impurity_radius, ss.max_expected_impurities_in_cell);
        QueryPerformanceCounter(&endClock);
        CreateMetricsLogs(i, GetElapsedTime(), grid);

        sample_results_coh[i] = RunSample(es, ss, i, true,  grid);
        sample_results_inc[i] = RunSample(es, ss, i, false, grid);

        std::cout << std::endl;
    }

    FinishResults(sample_results_coh, sample_results_inc);
}

SampleResult SimulationRunner::RunSample(Simulation& es, const Settings &settings, const int sample_index, const bool coherent, const Grid& grid)
{
    const int T = cfg.temperatures.size();
    const int N = cfg.magnetic_fields.size();
    SampleResult sr(T, N);

    auto metrics_path = GetMetricsPath(sample_index);
    double nlifetimes = pow(cfg.particles_per_row - 1, 2) * 4.0 * cfg.quadrant_integral_steps;

    es.InitSample(grid, settings, coherent);

    SampleMetrics sample_metrics(coherent, N, nlifetimes, grid.GetCellsPerRow(), grid.GetUniqueImpurityCount());
    
    for (int j = 0; j < N; j++) {
        Metrics metrics(j, nlifetimes, grid.GetCellsPerRow(), grid.GetUniqueImpurityCount());
    
        // Main compute method.
        QueryPerformanceCounter(&beginClock);
        es.ComputeLifetimes(cfg.magnetic_fields[j], grid, metrics);
        QueryPerformanceCounter(&endClock);
        
        metrics.time_elapsed_lifetimes = GetElapsedTime();
        sample_metrics.iteration_metrics[j] = metrics;

        for (int i = 0; i < T; i++) {
            auto iteration = es.DeriveTemperature(cfg.temperatures[i]);
            sr.results[i][j] = iteration.result;

            Logger::LogImages(GetImagePath(i, j, sample_index, coherent), cfg.particles_per_row - 1, iteration);
        }

        std::cout << ".";
    }

    Logger::LogSampleMetrics(metrics_path, sample_metrics);

    return sr;
}

// Finish results by averaging all samples and logging the result.
void SimulationRunner::FinishResults(const std::vector<SampleResult> sample_results_coh, const std::vector<SampleResult> sample_results_inc)
{
    const int S = cfg.num_samples;

    for (int i = 0; i < cfg.temperatures.size(); i++) {
        double temperature = cfg.temperatures[i];

        for (int j = 0; j < cfg.magnetic_fields.size(); j++) {
            Sigma coherent, incoherent;
            double dxx_squared = 0;

            for (int s = 0; s < S; s++) {
                // Sommeer <MF,T> voor alle samples.
                auto coh = sample_results_coh[s].results[i][j];
                auto inc = sample_results_inc[s].results[i][j];

                coherent += coh;
                incoherent += inc;

                dxx_squared += pow(coh.xx + inc.xx, 2);
            }

            coherent.xx /= S;
            coherent.xy /= S;
            incoherent.xx /= S;
            incoherent.xy /= S;

            DataRow row;
            row.temperature = temperature;
            row.magnetic_field = cfg.magnetic_fields[j];
            row.coherent = coherent;
            row.coherent = incoherent;

            double sxx_sq_exp = dxx_squared / S + 1e-15;
            double sxx_exp = (coherent.xx + incoherent.xx) / S;
            double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (S - 1));
            row.xxd = sxx_std / sxx_exp;

            Logger::LogResult(GetResultPath(i), row);
        }
    }
}

void SimulationRunner::ParseConfig(std::string file)
{
    std::unordered_map<std::string, std::string> values;

    std::filebuf fb;
    if (!fb.open(file, std::ios::in)) {
        std::cout << "Could not load config file.";
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

    base_output_directory = GetAvailableDirectory(values.at("output_directory"));

    cfg.num_samples             = atoi(values.at("num_samples").c_str());
    cfg.quadrant_integral_steps = atoi(values.at("integrand_steps").c_str());
    cfg.particles_per_row       = atoi(values.at("dimension").c_str());;

    {
        Range magnetic_field;
        magnetic_field.min       = atof(values.at("magnetic_field_min").c_str());
        magnetic_field.max       = atof(values.at("magnetic_field_max").c_str());
        magnetic_field.n         = atoi(values.at("magnetic_field_n").c_str());
        magnetic_field.step_size = (magnetic_field.max - magnetic_field.min) / (double)(magnetic_field.n - 1);

        cfg.magnetic_fields.resize(magnetic_field.n);
        for (int i = 0; i < magnetic_field.n; i++) {
            cfg.magnetic_fields[i] = magnetic_field.min + i * magnetic_field.step_size;
        }
    }

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

    cfg.settings.particle_speed                  = atof(values.at("particle_speed").c_str());
    cfg.settings.tau                             = atof(values.at("tau").c_str());
    cfg.settings.alpha                           = atof(values.at("alpha").c_str());
    cfg.settings.is_clockwise                    = atoi(values.at("clockwise").c_str()) > 0;

    cfg.settings.region_size                     = atof(values.at("region_size").c_str());
    cfg.settings.region_extends                  = atof(values.at("region_extends").c_str());
    cfg.settings.impurity_density                = atof(values.at("impurity_density").c_str());
    cfg.settings.impurity_radius                 = atof(values.at("impurity_radius").c_str());
    cfg.settings.max_expected_impurities_in_cell = atoi(values.at("max_expected_impurities_in_cell").c_str());
}

void SimulationRunner::PrintSimulationInfo() const
{
    std::cout << "Starting simulation (" << cfg.num_samples << "samples)\n";
    std::cout << "\tWriting to:" << base_output_directory << "\n\n";

    std::cout << "/" << std::string(cfg.magnetic_fields.size() * 2, '-') << "\\\n";
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

void SimulationRunner::CreateMetricsLogs(const int sample_index, const double elapsed_time, const Grid& grid) const
{
    GlobalMetrics gm;
    gm.particles_per_row                  = cfg.particles_per_row - 1;
    gm.phi_values                         = 4 * cfg.quadrant_integral_steps;
    gm.cells_per_row                      = grid.GetCellsPerRow();
    gm.unique_impurity_count              = grid.GetUniqueImpurityCount();
    gm.additional_impurities              = grid.GetTotalImpurityCount() - gm.unique_impurity_count;
    gm.grid_time_elapsed                  = elapsed_time; 
    gm.avg_impurities_in_cell             = grid.GetTotalImpurityCount() / (double)(gm.cells_per_row * gm.cells_per_row);
    gm.avg_impurities_in_cell_overlapping = gm.avg_impurities_in_cell - cfg.settings.max_expected_impurities_in_cell;

    Logger::CreateMetricsLog(GetMetricsPath(sample_index), gm);
}
