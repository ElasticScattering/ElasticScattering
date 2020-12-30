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

    std::cout << "Starting simulation (" << cfg.num_samples << " samples)\n";
    std::cout << "   Output directory: " << cfg.output_directory << "\n\n";
    std::cout << "/" << std::string(cfg.magnetic_fields.size() * 2, '-') << "\\\n";

    CreateOutputDirectories();

    std::vector<SampleResult> sample_results_coh(cfg.num_samples);
    std::vector<SampleResult> sample_results_inc(cfg.num_samples);

    std::random_device random_device;
    SimulationCPU es(cfg.particles_per_row-1, cfg.quadrant_integral_steps);

    const Settings& ss = cfg.settings;

    for (int i = 0; i < cfg.num_samples; i++) {
        std::cout << " ";
        
        QueryPerformanceCounter(&beginClock);
        auto grid = Grid(random_device(), ss.region_size, ss.region_extends, ss.impurity_density, ss.impurity_radius, ss.max_expected_impurities_in_cell);
        QueryPerformanceCounter(&endClock);
        CreateMetricsLogs(i, GetElapsedTime(), grid);

        sample_results_inc[i] = RunSample(es, ss, i, false, grid);
        sample_results_coh[i] = RunSample(es, ss, i, true,  grid);

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

            coherent.xx   /= S;
            coherent.xy   /= S;
            incoherent.xx /= S;
            incoherent.xy /= S;

            DataRow row(temperature, cfg.magnetic_fields[j], coherent, incoherent);

            double sxx_sq_exp = dxx_squared / S + 1e-15;
            double sxx_exp = (coherent.xx + incoherent.xx) / S;
            double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (S - 1));
            row.xxd = sxx_std / sxx_exp;

            Logger::LogResult(GetResultPath(i), row);
        }
    }
}

void SimulationRunner::CreateOutputDirectories() const
{
    if (!std::filesystem::exists(cfg.base_output_directory))
        std::filesystem::create_directory(cfg.base_output_directory);

    std::filesystem::create_directory(cfg.output_directory);

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
