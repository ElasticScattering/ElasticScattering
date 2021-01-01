#include "SimulationRunner.h"

#include "sim/Simulation.h"
#include "sim/Grid.h"
#include "SimulationResult.h"
#include "Logger.h"

#include <string>
#include <filesystem>

#include <random>

void SimulationRunner::Run(const InitParameters& init)
{
    QueryPerformanceFrequency(&clockFrequency);

    cfg = SimulationConfiguration::ParseFromeFile(init.config_file);
    std::cout << "Starting simulation (" << init.config_file << ")";
    if (cfg.output_type != OutputType::Nothing) {
        std::cout << " | Output: " << cfg.output_directory << "\n";
        CreateOutputDirectory();
    }
    std::cout << std::endl;

    std::vector<SampleResult> sample_results_coh(cfg.num_samples);
    std::vector<SampleResult> sample_results_inc(cfg.num_samples);

    std::random_device random_device;
    SimulationCPU es(cfg.particles_per_row-1, cfg.quadrant_phi_steps);

    const Settings& ss = cfg.settings;

    LARGE_INTEGER beginGridClock, endGridClock;
    LARGE_INTEGER beginTotalClock, endTotalClock;
    QueryPerformanceCounter(&beginTotalClock);

    for (int i = 0; i < cfg.num_samples; i++) {
        if (cfg.output_type == OutputType::All)
            CreateSampleOutputDirectory(i);

        QueryPerformanceCounter(&beginGridClock);
        auto seed = random_device();
        auto grid = Grid(seed, ss.region_size, ss.region_extends, ss.impurity_density, ss.impurity_radius, ss.max_expected_impurities_in_cell);
        QueryPerformanceCounter(&endGridClock);
        double grid_creation_time = GetElapsedTime(beginGridClock, endGridClock);

        if (i == 0 && cfg.output_type != OutputType::Nothing) {
            GlobalMetrics gm;
            gm.particles_per_row = cfg.particles_per_row - 1;
            gm.phi_values = 4 * cfg.quadrant_phi_steps;
            gm.cells_per_row = grid.GetCellsPerRow();
            gm.unique_impurity_count = grid.GetUniqueImpurityCount();
            gm.grid_creation_time = grid_creation_time;

            Logger::CreateSampleMetricsLog(GetMetricsPath(), gm);
        }

        sample_results_coh[i] = RunSample(es, ss, i, true,  grid);
        sample_results_inc[i] = RunSample(es, ss, i, false, grid);

        if (cfg.output_type == OutputType::All)
            Logger::LogSampleResults(GetSampleResultsPath(i), sample_results_coh[i], sample_results_inc[i]);
    }

    FinishResults(sample_results_coh, sample_results_inc);

    QueryPerformanceCounter(&endTotalClock);
    std::cout << std::endl << "Finished! (" << GetElapsedTime(beginTotalClock, endTotalClock) << "s)" << std::endl;
}

SampleResult SimulationRunner::RunSample(Simulation& es, const Settings &settings, const int sample_index, const bool coherent, const Grid& grid)
{
    const int T = cfg.temperatures.size();
    const int N = cfg.magnetic_fields.size();
    SampleResult sr(T, N);

    auto metrics_path = GetMetricsPath();
    double nlifetimes = pow(cfg.particles_per_row - 1, 2) * 4.0 * cfg.quadrant_phi_steps;

    es.InitSample(grid, settings, coherent);

    SampleMetrics sample_metrics(sample_index, coherent, N, nlifetimes);
    sample_metrics.total_indexed_impurities = grid.GetTotalImpurityCount();
    sample_metrics.impurity_count           = grid.GetUniqueImpurityCount();
    sample_metrics.total_cells              = pow(grid.GetCellsPerRow(), 2);
    sample_metrics.seed                     = grid.GetSeed();

    auto raw_sample_string = std::to_string(sample_index + 1);
    auto sample_string = std::string(cfg.digits_in_sample_num - raw_sample_string.length(), '0') + raw_sample_string + (coherent ? " C" : " I") + " [";

    std::cout << '\r' << sample_string << std::string(cfg.magnetic_fields.size(), '.') << "]";
    std::cout << '\r' << sample_string;

    LARGE_INTEGER beginLifetimesClock, endLifetimesClock;

    for (int j = 0; j < N; j++) {
        Metrics metrics;
    
        // Main compute method.
        QueryPerformanceCounter(&beginLifetimesClock);
        es.ComputeLifetimes(cfg.magnetic_fields[j], grid, metrics);
        QueryPerformanceCounter(&endLifetimesClock);
        
        metrics.time_elapsed_lifetimes = GetElapsedTime(beginLifetimesClock, endLifetimesClock);
        metrics.real_lifetimes = sample_metrics.total_lifetimes - metrics.particles_inside_impurity;
        sample_metrics.iteration_metrics[j] = metrics;

        for (int i = 0; i < T; i++) {
            if (cfg.output_type == OutputType::Nothing) 
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

    if (cfg.output_type != OutputType::Nothing)
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

            if (cfg.num_samples == 1) {
                row.xxd = 0;
            }
            else {
                double sxx_sq_exp = dxx_squared / S + 1e-15;
                double sxx_exp = (coherent.xx + incoherent.xx) / S;
                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(S - 1));
                row.xxd = sxx_std / sxx_exp;
            }

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

