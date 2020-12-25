
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Logger.h"

#include "scattering/escl/constants.h"

#include <algorithm>
#include <filesystem>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip> 
#include <ctime>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

void Logger::CreateResultLog(std::string file_path, const SimulationConfiguration& cfg, double temperature)
{
    std::ofstream file;
    file.open(file_path);
    
    ScatteringParameters sp = cfg.scattering_params;

    file << "# Elastic Scattering simulation results summary." << std::endl;
    file << "# Simulations " << cfg.magnetic_field_range.n << std::endl;
    file << "# Samples     " << cfg.num_samples << std::endl;
    file << "# Particles   " << sp.dim-1 << "x" << sp.dim-1 << std::endl;
    file << "# Phi steps   " << sp.integrand_steps << std::endl;

    file << std::scientific << std::setprecision(3);
    file << std::endl;
    file << "###########################" << std::endl;
    file << "# Scattering parameters   #" << std::endl;
    file << "#" << std::endl;
    file << "#\t" << "Temperature    " << temperature << std::endl;
    file << "#\t" << "Tau            " << sp.tau << std::endl;
    file << "#\t" << "Alpha          " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed " << sp.particle_speed << std::endl;
    file << "#\t" << "Clockwise      " << (sp.is_clockwise ? "True" : "False") << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Count          " << sp.impurity_count << std::endl;
    file << "#\t" << "Region size    " << sp.region_size << std::endl;
    file << "#\t" << "Region extends " << sp.region_extends << std::endl;
    file << "#\t" << "Density        " << sp.impurity_density << std::endl;
    file << "#\t" << "Radius         " << sp.impurity_radius << std::endl;

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
    file << "magnetic_field     sigma_xx_inc        sigma_xx_coh       sigma_xy_inc        sigma_xy_coh       delta_xx" << std::endl;
}

void Logger::CreateMetricsLog(const std::string file_path, const GlobalMetrics& gm)
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
    file << "Avg. imps/cell " << gm.avg_impurities_in_cell << " (" << gm.avg_impuritiies_in_cell_overlapping << " from overlapping)" << std::endl;
    file << std::fixed << std::setprecision(4);
    file << "Index time     " << gm.grid_time_elapsed << " seconds" << std::endl;
    file << std::endl << std::endl << std::endl << std::endl;
}

void Logger::LogMetrics(const std::string file_path, const Metrics& metrics)
{
    auto f = fopen(file_path.c_str(), "a");
    const int metric_width = 26;
    const int tabbed_width = metric_width - 4;

    std::string label = "Metric (MF: " + std::to_string(metrics.magnetic_field_index) + ")";

    fprintf(f, "%-*s| %-13s\n", metric_width, label.c_str(), "Value");
    fprintf(f, "%s|%s\n", std::string(metric_width, '-').c_str(), std::string(14, '-').c_str());
    fprintf(f, "%-*s|\n", metric_width, "Time spent");
    fprintf(f, "\t%-*s| %.3fs\n", tabbed_width, "Computing lifetimes", metrics.time_elapsed_lifetimes);
    fprintf(f, "\t%-*s| %.3fs\n", tabbed_width, "Rest", metrics.time_elapsed_temperatures);
    fprintf(f, "%-*s|\n", metric_width, "");
    fprintf(f, "%-*s| %f mil\n", metric_width, "Impurity intersections", metrics.mln_impurity_intersections);
    fprintf(f, "\t%-*s| %.2f\n", tabbed_width, "Particle average", metrics.prt_impurity_intersections);
    fprintf(f, "\t%-*s| %.2f%%\n", tabbed_width, "Particle percentage", metrics.pct_prt_impurity_intersections);
    fprintf(f, "%-*s|\n", metric_width, "");
    fprintf(f, "%-*s| %f mil\n", metric_width, "Cells passed", metrics.mln_cells_passed);
    fprintf(f, "\t%-*s| %.2f\n", tabbed_width, "Particle average", metrics.prt_cells_passed);
    fprintf(f, "\t%-*s| %.2f%%\n", tabbed_width, "Particle percentage", metrics.pct_prt_cells_passed);
    fprintf(f, "%-*s|\n", metric_width, "");
    fprintf(f, "%-*s| %i\n", metric_width, "Escaped", metrics.particles_escaped);
    fprintf(f, "%-*s| %.2f%%\n", metric_width, "Started inside impurity", metrics.pct_particles_inside_impurity);
    fprintf(f, "%s\n", std::string(metric_width+15, '-').c_str());
    fprintf(f, "\n\n");

    fclose(f);
}

void Logger::LogResult(const std::string file_path, const DataRow row)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    std::string s = "   ";
    file << std::scientific << std::setprecision(10);

    file << row.magnetic_field << s << row.incoherent.xx << s << row.coherent.xx << s << row.incoherent.xy << s << row.coherent.xy << s << row.xxd << std::endl;
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

void Logger::LogImages(const std::string file_path, const int dim, const IterationResult iteration)
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
