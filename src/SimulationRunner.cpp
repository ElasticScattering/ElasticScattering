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
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    ParseConfig(init.config_file);
    PrintSimulationInfo();

    CreateOutputDirectories();

    SimulationCPU es;

    for (int i = 0; i < cfg.magnetic_field_range.n; i++) {
        QueryPerformanceCounter(&beginClock);
        RunIteration(i, es);
        QueryPerformanceCounter(&endClock);

        double time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

        printf("Simulation completed (%d/%d) - %f\n", i + 1, cfg.magnetic_field_range.n, time_elapsed);
    }
}

#define GetElapsedTime(s) ((double)(s.endClock.QuadPart - s.beginClock.QuadPart) / s.clockFrequency.QuadPart)

#define RunSample(sp, total_results, sample_results, iteration_metrics)           \
std::vector<Sigma> sample_results(N);                                   \
{                                                                       \
    iteration_metrics.cells_passed = 0; \
    iteration_metrics.impurity_intersections = 0; \
    iteration_metrics.particles_inside_impurity = 0; \
    iteration_metrics.particles_escaped = 0; \
    QueryPerformanceCounter(&sample_metrics.beginClock);                \
    es.ComputeLifetimes(sp, grid, iteration_metrics);                             \
    QueryPerformanceCounter(&sample_metrics.endClock);                  \
                                                                        \
    iteration_metrics.time_elapsed_lifetimes = GetElapsedTime(sample_metrics);    \
                                                                        \
    QueryPerformanceCounter(&sample_metrics.beginClock);                \
    for (int i = 0; i < cfg.temperatures.size(); i++) {                 \
        UpdateTemperature(sp, cfg.temperatures[i]);                     \
        auto iteration = es.DeriveTemperature(cfg.temperatures[i]);     \
                                                                        \
        total_results[i] += iteration.result;                           \
        sample_results[i] = iteration.result;                           \
                                                                        \
        auto image_path = GetImagePath(i, magnetic_index, j, sp.is_incoherent == 1); \
        Logger::LogImages(image_path, sp.dim - 1, iteration);           \
    }                                                                   \
    QueryPerformanceCounter(&sample_metrics.endClock);                  \
    iteration_metrics.time_elapsed_temperatures = GetElapsedTime(sample_metrics); \
	                                                                    \
	double nparticles = pow(sp.dim-1, 2);                               \
    double n_lt = nparticles * 4.0 * sp.integrand_steps;                \
                                                                        \
    iteration_metrics.particles_inside_impurity /= (4.0 * sp.integrand_steps);    \
    iteration_metrics.pct_particles_inside_impurity = 100 * iteration_metrics.particles_inside_impurity / nparticles; \
                                                                                                  \
    iteration_metrics.pct_particles_escaped = 100 * iteration_metrics.particles_escaped / nparticles; \
                                                                                                                        \
    iteration_metrics.prt_cells_passed = (double)iteration_metrics.cells_passed / n_lt;                               \
    iteration_metrics.pct_prt_cells_passed = 100 * iteration_metrics.prt_cells_passed / pow(sp.cells_per_row, 2);     \
                                                                                                  \
    iteration_metrics.prt_impurity_intersections = (double)iteration_metrics.impurity_intersections / n_lt;           \
    iteration_metrics.pct_prt_impurity_intersections = 100 * iteration_metrics.prt_impurity_intersections / sp.impurity_count; \
}                                                                       

void SimulationRunner::RunIteration(const int magnetic_index, Simulation& es)
{
    UpdateMagneticField(cfg.scattering_params, cfg.magnetic_field_range.min + cfg.magnetic_field_range.step_size * magnetic_index);

    ScatteringParameters sp_inc = cfg.scattering_params;
    sp_inc.is_incoherent = 1;
    CompleteSimulationParameters(sp_inc);

    ScatteringParameters sp_coh = cfg.scattering_params;
    sp_coh.is_incoherent = 0;
    CompleteSimulationParameters(sp_coh);

    std::random_device random_device;

    const int N = cfg.temperatures.size();
    std::vector<Sigma>  coherent_results(N);
    std::vector<Sigma>  incoherent_results(N);
    std::vector<double> total_dxx_sq(N);

    Logger::CreateMetricsLog(GetMetricsPath(magnetic_index), cfg.scattering_params);

    for (int j = 0; j < cfg.samples_per_run; j++) {
        auto grid = Grid(sp_coh.impurity_count, random_device(), sp_coh.impurity_spawn_range, sp_coh.impurity_radius, sp_coh.cells_per_row);

        SampleMetrics sample_metrics;
        QueryPerformanceFrequency(&sample_metrics.clockFrequency);
        sample_metrics.actual_impurity_count = grid.GetImpurities().size();

        RunSample(sp_coh, coherent_results, sample_coherent_results, sample_metrics.coherent);
        RunSample(sp_inc, incoherent_results, sample_incoherent_results, sample_metrics.incoherent);

        Logger::LogMetrics(GetMetricsPath(magnetic_index), sample_metrics, j);

        for (int i = 0; i < N; i++) {
            total_dxx_sq[i] += pow(sample_coherent_results[i].xx + sample_incoherent_results[i].xx, 2);
        }

        printf("\tSample completed!\n");
    }

    double sample_count = (double)cfg.samples_per_run;
    for (int i = 0; i < N; i++)
    {
        coherent_results[i].xx /= sample_count;
        coherent_results[i].xy /= sample_count;

        incoherent_results[i].xx /= sample_count;
        incoherent_results[i].xy /= sample_count;
    }

    for (int i = 0; i < N; i++)
    {
        DataRow row;
        row.magnetic_field = cfg.scattering_params.magnetic_field;
        row.temperature = cfg.temperatures[i];
        row.coherent = coherent_results[i];
        row.incoherent = incoherent_results[i];

        double sxx_sq_exp = total_dxx_sq[i] / sample_count + 1e-15;
        double sxx_exp = (coherent_results[i].xx + incoherent_results[i].xx) / sample_count;
        double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (sample_count - 1));

        row.xxd = sxx_std / sxx_exp;

        Logger::LogResult(GetResultPath(i), row);
    }
}

void SimulationRunner::UpdateMagneticField(ScatteringParameters& sp, const double magnetic_field) {
    sp.magnetic_field = magnetic_field;
    sp.angular_speed = E * sp.magnetic_field / M;
}

void SimulationRunner::UpdateTemperature(ScatteringParameters& sp, const double temperature) {
    sp.temperature = temperature;

    if (sp.is_incoherent == 1) {
        sp.tau = HBAR / (KB * sp.temperature);
        sp.default_max_lifetime = 15.0 * sp.tau;
    }
}

void SimulationRunner::CompleteSimulationParameters(ScatteringParameters& sp) {
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

void SimulationRunner::PrintSimulationInfo() const
{
    printf("Starting simulation...\n");
    printf("\tWriting to: %s\n\n", cfg.base_output_directory.c_str());
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

    cfg.samples_per_run = atoi(values.at("samples_per_run").c_str());
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

    sp.is_clockwise = atoi(values.at("is_clockwise").c_str());

    CompleteSimulationParameters(sp);

    cfg.scattering_params = sp;
}

void SimulationRunner::CreateOutputDirectories() const
{
    std::filesystem::create_directory(cfg.base_output_directory + "/Metrics");

    for (int i = 0; i < cfg.temperatures.size(); i++) {
        Logger::CreateTemperatureLog(GetResultPath(i), cfg, cfg.temperatures[i]);

        auto temp_dir = cfg.base_output_directory + "/Images T" + std::to_string(i);
        std::filesystem::create_directory(temp_dir);

        for (int i = 0; i < cfg.magnetic_field_range.n; i++) {
            std::filesystem::create_directory(temp_dir + "/MF" + std::to_string(i));
        }
    }
}

std::string SimulationRunner::GetImagePathBase(int t_idx, int m_idx) const
{
    return cfg.base_output_directory + "/Images T" + std::to_string(t_idx) + "/MF" + std::to_string(m_idx);
}

std::string SimulationRunner::GetImagePath(int t_idx, int m_idx, int s_idx, bool incoherent) const
{
    auto type = incoherent ? " Incoherent" : " Coherent";
    return GetImagePathBase(t_idx, m_idx) + "/S" + std::to_string(s_idx) + type + ".png";
}

std::string SimulationRunner::GetMetricsPath(int m_idx) const
{
    return cfg.base_output_directory + "/Metrics/Metrics MF" + std::to_string(m_idx) + ".txt";
}

std::string SimulationRunner::GetResultPath(int t_idx) const
{
    return cfg.base_output_directory + "/T" + std::to_string(t_idx) + ".dat";
}
