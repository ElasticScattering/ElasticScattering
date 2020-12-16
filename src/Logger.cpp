
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

void Logger::CreateTemperatureLog(std::string file_path, const SimulationConfiguration& cfg, double temperature)
{
    std::ofstream file;
    file.open(file_path);
    
    ScatteringParameters sp = cfg.scattering_params;

    file << "# Elastic Scattering simulation results summary." << std::endl;
    file << "# Simulations " << cfg.magnetic_field_range.n << std::endl;
    file << "# Samples     " << cfg.samples_per_run << std::endl;
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
    file << "#\t" << "Clockwise      " << ((sp.is_clockwise == 1) ? "True" : "False") << std::endl;
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
    file << "#\n# Other:" << std::endl;
    file << "#" << "Grid cells       " << sp.cells_per_row << "x" << sp.cells_per_row << std::endl;
    file << "#" << "Grid cell size   " << sp.cell_size << std::endl;
    file << std::endl;
    file << "###########################" << std::endl;
    file << "# Results                 #" << std::endl;
    file << "#" << std::endl << std::endl;
    file << "magnetic_field     sigma_xx_inc        sigma_xx_coh       sigma_xy_inc        sigma_xy_coh       delta_xx" << std::endl;
}

void Logger::CreateMetricsLog(const std::string file_path, const ScatteringParameters& sp)
{
    std::ofstream file;
    file.open(file_path);

    int phi_values = 4 * sp.integrand_steps;
    int lifetimes = pow(sp.dim - 1, 2) * phi_values;
    
    file << "# Metrics collected during particle lifetime computation." << std::endl << std::endl;
    file << std::setprecision(4);
    file << ":/ Problem size" << std::endl;
    file << "Magnetic field " << sp.magnetic_field << std::endl;
    file << "Particles      " << sp.dim - 1 << "x" << sp.dim - 1 << std::endl;
    file << "Phi values     " << phi_values << std::endl;
    file << "Lifetimes      " << lifetimes << std::endl;
    file << "Impurities     " << sp.impurity_count << std::endl;
    file << "Cells          " << sp.cells_per_row << "x" << sp.cells_per_row << std::endl;
    file << std::endl << std::endl << std::endl;
}

void Logger::LogMetrics(const std::string file_path, const SampleMetrics& metrics, const int sample_id)
{
    /*
    Niet erg interessant maar wel per iteratie...
    double avg_impurities_in_cell = metrics.actual_impurity_count / (double)(sp.cells_per_row * sp.cells_per_row);
    file << "Avg. impurities in cell " << avg_impurities_in_cell << " (" << avg_impurities_in_cell - sp.max_expected_impurities_in_cell << "from overlapping)" << std::endl;
    */

    const int metric_width = 26;
    const int value_width = 15;

    auto f = fopen(file_path.c_str(), "a");
    fprintf(f, ":/ Sample %i\n", sample_id);
    fprintf(f, "\n");
    fprintf(f, "%-*s| %-13s| %s\n", metric_width, "Metric", "Coherent", "Incoherent");
    fprintf(f, "%s|%s|%s\n", std::string(metric_width, '-').c_str(), std::string(14, '-').c_str(), std::string(14, '-').c_str());
    fprintf(f, "%-*s| %-13f| %f\n", metric_width, "Time passed LT (s)", metrics.coherent.time_elapsed_lifetimes, metrics.incoherent.time_elapsed_lifetimes);
    fprintf(f, "%-*s| %-13f| %f\n", metric_width, "Time passed Rest (s)", metrics.coherent.time_elapsed_temperatures, metrics.incoherent.time_elapsed_temperatures);
    fprintf(f, "%-*s| %-13s|\n", metric_width, "", "");
    fprintf(f, "%-*s| %-13i| %i\n", metric_width, "Impurity intersections", metrics.coherent.impurity_intersections, metrics.incoherent.impurity_intersections);
    fprintf(f, "\t%-*s| %-13f| %f\n", metric_width - 4, "Particle average", metrics.coherent.prt_impurity_intersections, metrics.incoherent.prt_impurity_intersections);
    fprintf(f, "\t%-*s| %-13f| %f\n", metric_width - 4, "Particle percentage", metrics.coherent.pct_prt_impurity_intersections, metrics.incoherent.pct_prt_impurity_intersections);
    fprintf(f, "%-*s| %-13s|\n", metric_width, "", "");
    fprintf(f, "%-*s| %-13i| %i\n", metric_width, "Cells passed", metrics.coherent.cells_passed, metrics.incoherent.cells_passed);
    fprintf(f, "\t%-*s| %-13f| %f\n", metric_width - 4, "Particle average", metrics.coherent.prt_cells_passed, metrics.incoherent.prt_cells_passed);
    fprintf(f, "\t%-*s| %-13f| %f\n", metric_width - 4, "Particle percentage", metrics.coherent.pct_prt_cells_passed, metrics.incoherent.pct_prt_cells_passed);
    fprintf(f, "%-*s| %-13s|\n", metric_width, "", "");
    fprintf(f, "%-*s| %-13i| %i\n", metric_width, "Escaped", metrics.coherent.particles_escaped, metrics.incoherent.particles_escaped);
    fprintf(f, "%-*s| %-13i| %i\n", metric_width, "Started inside impurity", metrics.coherent.particles_inside_impurity, metrics.incoherent.particles_inside_impurity);
    fprintf(f, "\n");
    fprintf(f, "\n");

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

void Logger::FinishLog(const std::string file_path, const Metrics metrics)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    auto date_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    file << "\n\n\n###########################" << std::endl;
    file << "# Statistics              #" << std::endl;
    file << "#" << std::endl << std::endl;
    file << "# Completed on " << std::put_time(std::localtime(&date_time), "%F %T") << "." << std::endl;
    file << "# Elapsed time " << metrics.time_elapsed_lifetimes << " seconds." << std::endl;
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
