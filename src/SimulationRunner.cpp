#include "SimulationRunner.h"

#include "scattering/Simulation.h"
#include "scattering/Grid.h"
#include "SimulationResult.h"
#include "Logger.h"

#include <string>
#include <filesystem>

#include <random>

void SimulationRunner::Run(const InitParameters& init)
{
    QueryPerformanceFrequency(&clockFrequency);

    cfg = SimulationConfiguration::ParseFromeFile(init.config_file);
    std::cout << "Starting simulation (" << cfg.num_samples << " samples) | Output: " << cfg.output_directory << "\n\n";
    CreateOutputDirectory();

    std::vector<SampleResult> sample_results_coh(cfg.num_samples);
    std::vector<SampleResult> sample_results_inc(cfg.num_samples);

    std::random_device random_device;
    SimulationCPU es(cfg.particles_per_row-1, cfg.quadrant_integral_steps);

    const Settings& ss = cfg.settings;

    for (int i = 0; i < cfg.num_samples; i++) {
        CreateSampleOutputDirectory(i);


        QueryPerformanceCounter(&beginClock);
        auto seed = random_device();
        auto grid = Grid(seed, ss.region_size, ss.region_extends, ss.impurity_density, ss.impurity_radius, ss.max_expected_impurities_in_cell);
        QueryPerformanceCounter(&endClock);
        
        if (cfg.logging_level == LoggingLevel::Everything)
            CreateMetricsLogs(i, GetElapsedTime(), grid);

        sample_results_inc[i] = RunSample(es, ss, i, false, grid);
        sample_results_coh[i] = RunSample(es, ss, i, true,  grid);

        if (cfg.logging_level == LoggingLevel::Everything)
            Logger::LogSampleResults(GetSampleResultsPath(i), sample_results_coh[i], sample_results_inc[i]);
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

    auto raw_sample_string = std::to_string(sample_index + 1);
    auto sample_string = std::string(cfg.digits_in_sample_num - raw_sample_string.length(), '0') + raw_sample_string + (coherent ? " C" : " I") + " [";

    std::cout << '\r' << sample_string << std::string(cfg.magnetic_fields.size(), '.') << "]";
    std::cout << '\r' << sample_string;

    for (int j = 0; j < N; j++) {
        Metrics metrics(j, nlifetimes, grid.GetCellsPerRow(), grid.GetUniqueImpurityCount());
    
        // Main compute method.
        QueryPerformanceCounter(&beginClock);
        es.ComputeLifetimes(cfg.magnetic_fields[j], grid, metrics);
        QueryPerformanceCounter(&endClock);
        
        metrics.time_elapsed_lifetimes = GetElapsedTime();
        sample_metrics.iteration_metrics[j] = metrics;

        for (int i = 0; i < T; i++) {
            if (cfg.logging_level == LoggingLevel::Silent) 
            {
                sr.results[i][j] = es.DeriveTemperature(cfg.temperatures[i]);
            }
            else 
            {
                auto iteration = es.DeriveTemperatureWithImages(cfg.temperatures[i]);
                sr.results[i][j] = iteration.result;

                Logger::LogImages(GetImagePath(i, j, sample_index, coherent), cfg.particles_per_row - 1, iteration);
            }
        }

        std::cout << "x";
    }

    if (cfg.logging_level == LoggingLevel::Everything)
        Logger::LogSampleMetrics(metrics_path, sample_metrics);

    return sr;
}

void SimulationRunner::FinishResults(const std::vector<SampleResult> sample_results_coh, const std::vector<SampleResult> sample_results_inc)
{
    const double S = (double)cfg.num_samples;

    for (int i = 0; i < cfg.temperatures.size(); i++) {
        double temperature = cfg.temperatures[i];
        Logger::CreateResultLog(GetResultPath(i), cfg, temperature);

        for (int j = 0; j < cfg.magnetic_fields.size(); j++) {
            Sigma coherent, incoherent;
            double dxx_squared = 0;

            for (int s = 0; s < cfg.num_samples; s++) {
                // Sommeer <MF,T> voor alle samples.
                auto coh = sample_results_coh[s].results[i][j];
                auto inc = sample_results_inc[s].results[i][j];

                coherent += coh;
                incoherent += inc;

                dxx_squared += pow(coh.xx + inc.xx, 2);
            }

            coherent.xx   /= S;
            coherent.xy   /= S;
            incoherent.xx /= S;
            incoherent.xy /= S;

            DataRow row(temperature, cfg.magnetic_fields[j], coherent, incoherent);

            double sxx_sq_exp = dxx_squared / S + 1e-15;
            double sxx_exp = (coherent.xx + incoherent.xx) / S;
            double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(S - 1));
            row.xxd = sxx_std / sxx_exp;

            Logger::LogResult(GetResultPath(i), row);
        }
    }
}

void SimulationRunner::CreateSampleOutputDirectory(const int sample_index) const
{
    auto sample_path = GetSamplePath(sample_index);
    std::filesystem::create_directory(sample_path);
    std::filesystem::create_directory(sample_path + "/Incoherent/");
    std::filesystem::create_directory(sample_path + "/Coherent/");

    Logger::CreateSampleResultLog(GetSampleResultsPath(sample_index), cfg);
}

void SimulationRunner::CreateOutputDirectory() const
{

    if (!std::filesystem::exists(cfg.base_output_directory))
        std::filesystem::create_directory(cfg.base_output_directory);

    std::filesystem::create_directory(cfg.output_directory);
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
