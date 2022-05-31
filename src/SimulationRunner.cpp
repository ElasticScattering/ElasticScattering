/* -------------------------------------------------------------------------
    This code is part of ElasticScattering.

    Copyright(C) 2022 Stijn Hinlopen

    This program is free software : you can redistribute it and /or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

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
    std::cout << "Starting simulation (" << init.config_file << ")" << std::endl;
    if (cfg.output_type != OutputType::Nothing) 
    {
        std::cout << "Output directory: \x1B[95m" << cfg.output_directory << "\x1B[0m" << std::endl;
        CreateOutputDirectory();
    }
    std::cout << std::endl;

    LARGE_INTEGER beginTotalClock, endTotalClock;
    QueryPerformanceCounter(&beginTotalClock);
    auto sr = RunSimulation();
    QueryPerformanceCounter(&endTotalClock);
    std::cout << std::endl << "Finished! (" << GetElapsedTime(beginTotalClock, endTotalClock) << "s)" << std::endl;

    FinishResults(sr);
}

SimulationResult SimulationRunner::RunSimulation() const
{
    const Settings& ss = cfg.settings;

    LARGE_INTEGER beginGridClock, endGridClock;
    QueryPerformanceCounter(&beginGridClock);
    auto grid = Grid(cfg.sample_seeds[0], ss.region_size, ss.region_extends, ss.impurity_density, ss.impurity_radius, ss.target_cell_population);
    QueryPerformanceCounter(&endGridClock);
    double grid_creation_time = GetElapsedTime(beginGridClock, endGridClock);

    if (cfg.output_type != OutputType::Nothing) {
        GlobalMetrics gm;
        gm.particles_per_row     = cfg.positions_per_row;
        gm.phi_steps             = 4 * cfg.particles_per_quadrant;
        gm.cells_per_row         = grid.GetCellsPerRow();
        gm.unique_impurity_count = grid.GetUniqueImpurityCount();
        gm.grid_creation_time    = grid_creation_time;

        Logger::CreateSampleMetricsLog(GetMetricsPath(), gm);
    }

    GridInformation gi;
    gi.cell_size              = grid.GetSettings().cell_size;
    gi.indexed_impurity_count = grid.GetTotalImpurityCount();
    gi.index_size             = grid.GetIndex().size();
    gi.region_size            = ss.region_size;
    gi.region_extends         = ss.region_extends;

    Simulation* es;
    if (cfg.use_gpu) es = new SimulationCL(cfg.positions_per_row, cfg.particles_per_quadrant, gi, cfg.print_info);
    else             es = new SimulationCPU(cfg.positions_per_row, cfg.particles_per_quadrant, gi);

    SimulationResult sr(cfg.num_samples);

    for (int i = 0; i < cfg.num_samples; i++) {
        if (cfg.output_type == OutputType::All) CreateSampleOutputDirectory(i);

        sr.coherent.push_back(  RunSample(*es, ss, i, true,  grid));
        sr.incoherent.push_back(RunSample(*es, ss, i, false, grid));

        if (cfg.output_type == OutputType::All) Logger::LogSampleResults(GetSampleResultsPath(i), sr.coherent[i], sr.incoherent[i]);

        if ( i < (cfg.num_samples-1))
            grid = Grid(cfg.sample_seeds[i+1], ss.region_size, ss.region_extends, ss.impurity_density, ss.impurity_radius, ss.target_cell_population);
    }

    return sr;
}

SampleResult SimulationRunner::RunSample(Simulation& es, const Settings &settings, const int sample_index, const bool coherent, const Grid& grid) const
{
    const int T = cfg.temperatures.size();
    const int N = cfg.magnetic_fields.size();
    SampleResult sr(T, N);

    auto metrics_path = GetMetricsPath();
    double nlifetimes = pow(cfg.positions_per_row, 2) * 4.0 * cfg.particles_per_quadrant;

    es.InitSample(grid, settings, coherent);

    SampleMetrics sample_metrics(sample_index, coherent, nlifetimes);
    sample_metrics.total_indexed_impurities = grid.GetTotalImpurityCount();
    sample_metrics.impurity_count           = grid.GetUniqueImpurityCount();
    sample_metrics.total_cells              = pow(grid.GetCellsPerRow(), 2);
    sample_metrics.seed                     = grid.GetSeed();
    sample_metrics.iteration_metrics.reserve(N);

    {
        auto raw_sample_string = std::to_string(sample_index + 1);
        auto sample_string = std::string(cfg.digits_in_sample_num - raw_sample_string.length(), '0') + raw_sample_string + (coherent ? " C" : " I") + " [";

        std::cout << '\r' << sample_string << std::string(cfg.magnetic_fields.size(), '.') << "]";
        std::cout << '\r' << sample_string;
    }

    for (int i = 0; i < N; i++) {
        if (cfg.output_type != OutputType::All)
        {
            sr.results[i] = es.ComputeSigmas(cfg.magnetic_fields[i], cfg.temperatures, grid, sample_metrics);
        }
        else
        {
            auto iteration = es.ComputeSigmasWithImages(cfg.magnetic_fields[i], cfg.temperatures, grid, sample_metrics);

            for (int j = 0; j < cfg.temperatures.size(); j++) {
                sr.results[i][j] = iteration[j].result;
                Logger::LogImages(GetImagePath(j, i, sample_index, coherent), cfg.positions_per_row, iteration[j]);
            }
        }

        std::cout << "x";
    }

    if (cfg.output_type != OutputType::Nothing)
        Logger::LogSampleMetrics(metrics_path, sample_metrics);

    return sr;
}

void SimulationRunner::FinishResults(const SimulationResult& sr) const
{
    const double S = (double)cfg.num_samples;

    for (int j = 0; j < cfg.temperatures.size(); j++) {
        double temperature = cfg.temperatures[j];

        auto current_file = GetResultPath(j);
        Logger::CreateResultLog(current_file, cfg, temperature);

        for (int i = 0; i < cfg.magnetic_fields.size(); i++) {
            Sigma coherent, incoherent;
            double dxx_squared = 0;

            for (int s = 0; s < cfg.num_samples; s++) {
                auto coh = sr.coherent  [s].results[i][j];
                auto inc = sr.incoherent[s].results[i][j];

                coherent    += coh;
                incoherent  += inc;
                dxx_squared += pow(coh.xx + inc.xx, 2);
            }

            coherent.xx   /= S;
            coherent.xy   /= S;
            incoherent.xx /= S;
            incoherent.xy /= S;

            DataRow row(temperature, cfg.magnetic_fields[i], coherent, incoherent);

            if (cfg.num_samples == 1) 
            {
                row.xxd = 0;
            }
            else 
            {
                double sxx_sq_exp = dxx_squared / S;
                double sxx_exp = coherent.xx + incoherent.xx;
                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (S - 1.0));
                row.xxd = sxx_std / sxx_exp;
            }

            Logger::LogResult(current_file, row);
        }
    }
}

void SimulationRunner::CreateSampleOutputDirectory(const int sample_index) const
{
    auto sample_path = GetSamplePath(sample_index);
    std::filesystem::create_directory(sample_path);
    std::filesystem::create_directory(sample_path + "/Incoherent/");
    std::filesystem::create_directory(sample_path + "/Coherent/");

    Logger::CreateSampleResultLog(GetSampleResultsPath(sample_index));
}

void SimulationRunner::CreateOutputDirectory() const
{
    if (!std::filesystem::exists(cfg.base_output_directory))
        std::filesystem::create_directory(cfg.base_output_directory);

    std::filesystem::create_directory(cfg.output_directory);

    auto slash_pos = cfg.config_path.rfind('/');
    auto config_name = slash_pos != std::string::npos ? cfg.config_path.substr(slash_pos, cfg.config_path.size() - slash_pos) : cfg.config_path;

    std::filesystem::copy_file(cfg.config_path, cfg.output_directory + "/" + config_name);
}
