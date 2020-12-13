
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

void Logger::CreateLog(const std::string result_file, const SimulationConfiguration& cfg, double temperature)
{
    std::ofstream file;
    file.open(result_file);

    file << "# Elastic Scattering simulation results" << std::endl;
    file << "# Number of runs: " << cfg.magnetic_field_range.n << std::endl;
    file << "# Samples per row: " << cfg.samples_per_run << std::endl;

    ScatteringParameters sp = cfg.scattering_params;

    file << std::scientific << std::setprecision(3);
    file << "#" << std::endl;
    file << "###########################" << std::endl;
    file << "# Scattering parameters:  #" << std::endl;
    file << "###########################" << std::endl;
    file << "#\t" << "Dimension:        " << sp.dim << std::endl;
    file << "#\t" << "Integrand steps:  " << sp.integrand_steps << std::endl;
    file << "#" << std::endl;
    file << "#\t" << "Temperature:      " << temperature << std::endl;
    file << "#\t" << "Tau:              " << sp.tau << std::endl;
    file << "#\t" << "Alpha:            " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed:   " << sp.particle_speed << std::endl;
    file << "#\t" << "Clockwise:        " << ((sp.is_clockwise == 1) ? "True" : "False") << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Count:            " << sp.impurity_count << std::endl;
    file << "#\t" << "Region size:      " << sp.region_size << std::endl;
    file << "#\t" << "Region extends:   " << sp.region_extends << std::endl;
    file << "#\t" << "Density:          " << sp.impurity_density << std::endl;
    file << "#\t" << "Radius:           " << sp.impurity_radius << std::endl;

    file << "#\n# Constants:" << std::endl;
    file << "#\t" << "Particle mass:    " << M << std::endl;
    file << "#\t" << "E:                " << E << std::endl;
    file << "#\t" << "HBAR:             " << HBAR << std::endl;
    file << "#\t" << "C:                " << C1 << std::endl;
    file << "#\t" << "KB:               " << KB << std::endl;
    file << "#" << std::endl;
    file << "###########################" << std::endl;
    file << "# Results:                #" << std::endl;
    file << "###########################" << std::endl << std::endl;
    file << "magnetic_field sigma_xx_inc sigma_xx_coh sigma_xy_inc sigma_xy_coh delta_xx" << std::endl;
}

void Logger::LogResult(const std::string file_path, const DataRow row)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    file << std::scientific << std::setprecision(10);
    file << row.magnetic_field << " " << row.incoherent.xx << " " << row.coherent.xx << " " << row.incoherent.xy << " " << row.coherent.xy << " " << row.xxd << std::endl;
}

void Logger::FinishLog(const std::string file_path, const double time_elapsed)
{
    std::ofstream file;
    file.open(file_path, std::ios_base::app);

    auto date_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    file << "\n\n\n###################################" << std::endl;
    file << "# Completed on: " << std::put_time(std::localtime(&date_time), "%F %T") << "." << std::endl;
    file << "# Elapsed time: " << time_elapsed << " seconds." << std::endl;
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

void Logger::LogImages(const std::string file, const int dim, const double scale, const IterationResult iteration)
{
    const int image_width  = dim * 3;
    const int image_height = dim;
    const int num_channels = 3;

    std::vector<unsigned char> pixels(image_width * image_height * num_channels);

    WriteImageSection(pixels, iteration.particle_lifetimes, dim, 0);
    WriteImageSection(pixels, iteration.sigmas.xx_buffer,   dim, 1);
    WriteImageSection(pixels, iteration.sigmas.xy_buffer,   dim, 2, true);

    stbi_write_png(file.c_str(), image_width, image_height, 3, pixels.data(), image_width * num_channels);
}
