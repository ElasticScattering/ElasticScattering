﻿
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Logger.h"

#include "scattering/escl/constants.h"

#include <algorithm>
#include <filesystem>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iomanip> 
#include <ctime>

#include <clocale>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

void Logger::CreateResultLog(const std::string file_path, const SimulationConfiguration& cfg, double temperature)
{
    std::ofstream file;
    file.open(file_path);

    const auto ss = cfg.settings;
    
    file << "# Elastic Scattering simulation results summary." << std::endl;
    file << "# Simulations " << cfg.magnetic_fields.size() << std::endl;
    file << "# Samples     " << cfg.num_samples << std::endl;
    file << "# Particles   " << cfg.particles_per_row-1 << "x" << cfg.particles_per_row-1 << std::endl;
    file << "# Phi steps   " << cfg.quadrant_integral_steps << std::endl;

    file << std::scientific << std::setprecision(3);
    file << std::endl;
    file << "###########################" << std::endl;
    file << "# Scattering parameters   #" << std::endl;
    file << "#" << std::endl;
    file << "#\t" << "Temperature    " << temperature << std::endl;
    file << "#\t" << "Tau            " << ss.tau << std::endl;
    file << "#\t" << "Alpha          " << ss.alpha << std::endl;
    file << "#\t" << "Particle speed " << ss.particle_speed << std::endl;
    file << "#\t" << "Clockwise      " << (ss.is_clockwise ? "True" : "False") << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Region size    " << ss.region_size << std::endl;
    file << "#\t" << "Region extends " << ss.region_extends << std::endl;
    file << "#\t" << "Density        " << ss.impurity_density << std::endl;
    file << "#\t" << "Radius         " << ss.impurity_radius << std::endl;

    file << "#\n# Constants:" << std::endl;
    file << "#\t" << "Particle mass  " << M << std::endl;
    file << "#\t" << "E              " << E << std::endl;
    file << "#\t" << "HBAR           " << HBAR << std::endl;
    file << "#\t" << "C              " << C1 << std::endl;
    file << "#\t" << "KB             " << KB << std::endl;
    file << std::endl;
    file << "###########################" << std::endl;
    file << "# Results                 #" << std::endl;
    file << "#" << std::endl << std::endl;
    file << "magnetic_field  sigma_xx_inc    sigma_xx_coh    sigma_xy_inc    sigma_xy_coh    delta_xx" << std::endl;
}


void Logger::LogResult(const std::string file_path, const DataRow& row)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    int p = 8;
    int w = p + 8;
    file << std::scientific << std::setprecision(p);

    file << std::setw(w) << std::left << row.magnetic_field 
         << std::setw(w) << std::left << row.incoherent.xx 
         << std::setw(w) << std::left << row.coherent.xx 
         << std::setw(w) << std::left << row.incoherent.xy
         << std::setw(w) << std::left << row.coherent.xy 
         << std::setw(w) << std::left << row.xxd << std::endl;
    
    //std::string s = "   ";
    //file << row.magnetic_field << s << row.incoherent.xx << s << row.coherent.xx << s << row.incoherent.xy << s << row.coherent.xy << s << row.xxd << std::endl;
}




void Logger::CreateSampleResultLog(const std::string file_path, const SimulationConfiguration& cfg)
{
    std::ofstream file;
    file.open(file_path);

    const auto ss = cfg.settings;

    file << "# Elastic Scattering sample results summary." << std::endl;
    file << "# Simulations " << cfg.magnetic_fields.size() << std::endl;
    file << "# Particles   " << cfg.particles_per_row - 1 << "x" << cfg.particles_per_row - 1 << std::endl;
    file << "# Phi steps   " << cfg.quadrant_integral_steps << std::endl;

    file << std::scientific << std::setprecision(3);
    file << std::endl;
    file << "###########################" << std::endl;
    file << "# Scattering parameters   #" << std::endl;
    file << "#" << std::endl;
    file << "#\t" << "Tau            " << ss.tau << std::endl;
    file << "#\t" << "Alpha          " << ss.alpha << std::endl;
    file << "#\t" << "Particle speed " << ss.particle_speed << std::endl;
    file << "#\t" << "Clockwise      " << (ss.is_clockwise ? "True" : "False") << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Region size    " << ss.region_size << std::endl;
    file << "#\t" << "Region extends " << ss.region_extends << std::endl;
    file << "#\t" << "Density        " << ss.impurity_density << std::endl;
    file << "#\t" << "Radius         " << ss.impurity_radius << std::endl;

    file << "###########################" << std::endl;
    file << "# Sample Results          #" << std::endl;
    file << "#" << std::endl << std::endl;

    int w = 16;
    file << std::setw(4) << std::left << "T"
         << std::setw(4) << std::left << "B"
         << std::setw(w) << std::left << "sigma_xx_inc"
         << std::setw(w) << std::left << "sigma_xx_coh"
         << std::setw(w) << std::left << "sigma_xy_inc" 
         << std::setw(w) << std::left << "sigma_xy_coh" << std::endl;
}

void Logger::LogSampleResults(const std::string file_path, const SampleResult coherent, const SampleResult incoherent)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    int p = 8;
    int w = p + 8;
    file << std::scientific << std::setprecision(p);

    for (int j = 0; j < coherent.results.size(); j++) {
        auto coh = coherent.results[j];
        auto inc = incoherent.results[j];

        for (int i = 0; i < coherent.results.size(); i++) {
            file << std::setw(4) << std::left << j
                 << std::setw(4) << std::left << i
                 << std::setw(w) << std::left << inc[i].xx
                 << std::setw(w) << std::left << coh[i].xx
                 << std::setw(w) << std::left << inc[i].xy
                 << std::setw(w) << std::left << coh[i].xy << std::endl;
        }
        file << std::endl;
    }
}







void Logger::CreateSampleMetricsLog(const std::string file_path, const GlobalMetrics& gm)
{
    std::ofstream file;
    file.open(file_path);

    int lifetimes = (int)pow(gm.particles_per_row, 2) * gm.phi_values;
    
    file << "# Metrics collected during particle lifetime computation." << std::endl << std::endl;
    file << ":/ Problem size" << std::endl;
    file << "Particles      " << gm.particles_per_row << " x " << gm.particles_per_row << std::endl;
    file << "Phi values     " << gm.phi_values << std::endl;
    file << "Lifetimes      " << lifetimes << std::endl;
    file << "Impurities     " << gm.unique_impurity_count << " (+ " << gm.additional_impurities << " from overlap)" << std::endl;
    file << "Cells          " << gm.cells_per_row << " x " << gm.cells_per_row << std::endl;
    file << "Avg. imps/cell " << gm.avg_impurities_in_cell << " (" << gm.avg_impurities_in_cell_overlapping << " from overlapping)" << std::endl;
    file << "Impurity seed  " << gm.seed << std::endl;
    file << std::fixed << std::setprecision(4);
    file << "Index time     " << gm.grid_time_elapsed * 1000.0 << "ms" << std::endl;
    file << std::endl;
}

void Logger::LogSampleMetrics(const std::string file_path, const SampleMetrics& sample_metrics)
{
    std::wofstream file;
    file.open(file_path, std::ios_base::app);

    const int metric_width = 26;
    const int value_width = 12;
    const auto metrics = sample_metrics.iteration_metrics;

    const auto loc = std::locale("en_US.UTF-8");
    file.imbue(loc);
    //const std::locale utf8_locale = std::locale(std::locale(), new std::codecvt_utf8<wchar_t>());
    
    file << std::endl << std::endl;

    // Header
    file << L'┌' << std::wstring(metric_width, L'─');
    for (int i = 0; i < metrics.size(); i++) file << L'┬' << std::wstring(value_width, L'─');
    file << L'┐' << std::endl;

    auto label = (sample_metrics.coherent) ? " Metric (Coherent)" : " Metric (Incoherent)";

    file << L'│' << std::setw(metric_width) << std::left << label;
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width) << std::left << std::wstring(L" MF ") + std::to_wstring(i);
    file << L'│' << std::endl;

    file << L'├' << std::wstring(metric_width, L'─');
    for (int i = 0; i < metrics.size(); i++) file << L'┼' << std::wstring(value_width, L'─');
    file << L'┤' << std::endl;

    // Rows
    file << L'│' << std::setw(metric_width) << std::left << " Avg. particle lifetime";
    for (int i = 0; i < metrics.size(); i++) file << L'│'  << std::setw(value_width) << std::right << std::setprecision(3) << metrics[i].avg_particle_lifetime;
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width) << std::left << " Time spent on lifetimes";
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width - 2) << std::right << std::setprecision(3) << (metrics[i].time_elapsed_lifetimes*1000.0) << "ms";
    file << L'│' << std::endl;

    file << L'│' << std::wstring(metric_width, ' '); // Empty line.
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::wstring(value_width, ' ');
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width) << std::left << " Total impurities tested";
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width - 4) << std::right << std::setprecision(3) << ((double)metrics[i].impurities_tested / 1'000'000.0) << " mil";
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width - 2) << std::left << "\t Average";
    for (int i = 0; i < metrics.size(); i++)
    {
        file << L'│' << std::setw(value_width) << std::right << std::setprecision(2) 
            << ((double)metrics[i].impurities_tested / (double)metrics[i].real_lifetimes);
    }
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width - 2) << std::left << "\t Average %";
    for (int i = 0; i < metrics.size(); i++)
    {
        file << L'│' << std::setw(value_width-1) << std::right << std::setprecision(2)
            << (100 * ((double)metrics[i].impurities_tested / (double)metrics[i].real_lifetimes) / (double)sample_metrics.impurity_count) << "%";
    }
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width - 2) << std::left << "\t Max";
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width) << std::right << std::setprecision(4) <<  (double)metrics[i].max_impurities_tested;
    file << L'│' << std::endl;

    file << L'│' << std::wstring(metric_width, ' '); // Empty line.
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::wstring(value_width, ' ');
    file << L'│' << std::endl;

    // # Cells Passed.
    file << L'│' << std::setw(metric_width) << std::left << " Grid cells passed";
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width - 4) << std::right << std::setprecision(3) << ((double)metrics[i].cells_passed / 1'000'000.0) << " mil";
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width - 2) << std::left << "\t Average";
    for (int i = 0; i < metrics.size(); i++)
    {
        file << L'│' << std::setw(value_width) << std::right << std::setprecision(2)
            << ((double)metrics[i].cells_passed / (double)metrics[i].real_lifetimes);
    }
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width - 2) << std::left << "\t Average %";
    for (int i = 0; i < metrics.size(); i++)
    {
        file << L'│' << std::setw(value_width - 1) << std::right << std::setprecision(2)
            << (100 * ((double)metrics[i].cells_passed / (double)metrics[i].real_lifetimes) / pow(sample_metrics.cells_per_row, 2)) << "%";
    }
    file << L'│' << std::endl;

    file << L'│' << std::setw(metric_width - 2) << std::left << "\t Max";
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width) << std::right << (double)metrics[i].max_cells_passed;
    file << L'│' << std::endl;

    file << L'│' << std::wstring(metric_width, ' '); // Empty line.
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::wstring(value_width, ' ');
    file << L'│' << std::endl;

    // # Particles escaped.
    file << L'│' << std::setw(metric_width) << std::left << " Particles escaped";
    for (int i = 0; i < metrics.size(); i++) file << L'│' << std::setw(value_width) << std::right << (double)metrics[i].particles_escaped;
    file << L'│' << std::endl;

    // # Particles inside impurity.
    file << L'│' << std::setw(metric_width) << std::left << " Started inside impurity";
    for (int i = 0; i < metrics.size(); i++)
    {
        file << L'│' << std::setw(value_width - 1) << std::right << std::setprecision(2)
            << (100 * ((double)metrics[i].particles_inside_impurity / sample_metrics.nlifetimes)) << "%";
    }
    file << L'│' << std::endl;
    
    // Table end.
    file << L'└' << std::wstring(metric_width, L'─');
    for (int i = 0; i < metrics.size(); i++) file << L'┴' << std::wstring(value_width, L'─');
    file << L'┘' << std::endl;
}





void Logger::WriteImageSection(std::vector<unsigned char> &pixels, const std::vector<double> &values, const int dim, const int image_id, const bool colored)
{
    const int image_width = dim * 3 * 3;
    const int image_height = dim;
    const int offset = image_id * dim * 3;

    double min_element =  INF;
    double max_element = -INF;

    for (int i = 0; i < values.size(); i++) {
        auto k = values[i];
        if (k < min_element) min_element = k;
        if (k > max_element) max_element = k;
    }

    int i, j, idx = 0;
    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
        {
            idx = j * image_width + i * 3 + offset;
            
            if (!colored) {
                auto k = (unsigned char)((values[j * dim + i] - min_element) / max_element * 255.0);
                pixels[idx]     = k;
                pixels[idx + 1] = k;
                pixels[idx + 2] = k;
            }
            else {
                auto c = values[j * dim + i];
                if (c < 0) {
                    pixels[idx]     = 0;
                    pixels[idx + 1] = 0;
                    pixels[idx + 2] = (unsigned char)(c / min_element * 255.0);
                }
                else {
                    pixels[idx]     = (unsigned char)(c / max_element * 255.0);
                    pixels[idx + 1] = 0;
                    pixels[idx + 2] = 0;
                }
            }
        }
}

void Logger::LogImages(const std::string file_path, const int dim, const IterationResult &iteration)
{
    const int image_width  = dim * 3;
    const int image_height = dim;
    const int num_channels = 3;

    std::vector<unsigned char> pixels(image_width * image_height * num_channels);

    WriteImageSection(pixels, iteration.particle_lifetimes, dim, 0);
    WriteImageSection(pixels, iteration.sigmas.xx_buffer,   dim, 1);
    WriteImageSection(pixels, iteration.sigmas.xy_buffer,   dim, 2, true);

    stbi_write_png(file_path.c_str(), image_width, image_height, 3, pixels.data(), image_width * num_channels);
}