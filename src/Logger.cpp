
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Logger.h"

#include "scattering/escl/constants.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <string>
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
    file << "###########################" << std::endl;
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

void Logger::LogImages(const std::string file, const int dim, const double scale, const IterationResult iteration)
{
    const int image_width = dim * 3;
    const int image_height = dim;

    std::vector<unsigned char> pixels(image_width * image_height);

    int i, j = 0;
    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + i] = (unsigned char)(iteration.particle_lifetimes[j * dim + i] / scale * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + (i + dim)] = (unsigned char)(iteration.sigmas.xx_buffer[j * dim + i] / 3.0 * 255.0);

    for (j = 0; j < dim; j++)
        for (i = 0; i < dim; i++)
            pixels[j * image_width + (i + dim * 2)] = (unsigned char)(iteration.sigmas.xy_buffer[j * dim + i] / 1.5 * 255.0);

    stbi_write_png(file.c_str(), image_width, image_height, 1, pixels.data(), image_width);
}
