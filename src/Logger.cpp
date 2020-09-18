#include "src/Logger.h"

#include "src/SimulationParameters.h"

#include <iostream>
#include <fstream>

void Logger::LogArrays(const SimulationParameters& sp, const std::vector<double>& xs, const std::vector<double>& yss) const
{
    std::ofstream file;
    file.open("test1.dat");
    file << "# Simulation parameters:" << std::endl;
    file << "#\t" << "Mode:             " << sp.mode << std::endl;
    file << "#\t" << "Integrand steps:  " << sp.integrand_steps << std::endl;
    file << "#\t" << "Dimension:        " << sp.dim << std::endl;
    file << "#\t" << "Temperature:      " << sp.temperature << std::endl;
    file << "#\t" << "Tau:              " << sp.tau << std::endl;
    file << "#\t" << "Magnetic field:   " << sp.magnetic_field << std::endl;
    file << "#\t" << "Phi:              " << sp.phi << std::endl;
    file << "#\t" << "Alpha:            " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed:   " << sp.particle_speed << std::endl;
    file << "#\t" << "Angular speed:    " << sp.angular_speed << std::endl;
    file << "#\t" << "Impurity density: " << sp.impurity_density << std::endl;
    file << "#\t" << "Impurity count:   " << sp.impurity_count << std::endl;
    file << "#\t" << "Impurity radius:  " << sp.impurity_radius << std::endl;
    file << "#\t" << "Impurity seed:    " << sp.impurity_seed << std::endl;
    file << "#\t" << "Region size:      " << sp.region_size << std::endl;
    file << "#\t" << "Region extends:   " << sp.region_extends << std::endl;
    file << "#\t" << "Diag. regions:    " << sp.is_diag_regions << std::endl;
    file << "#\t" << "Clockwise:        " << sp.is_clockwise << std::endl;
    file << "#\t" << "Incoherent:       " << sp.is_incoherent << std::endl;

    file << "\n" << std::endl;
    file << "# Results:" << std::endl;

    file << "# [Magnetic field] [S-XX] [S-XY] [rho]" << std::endl;

    int n = xs.size();
    int stride = yss.size() / n;

    std::cout.precision(17);

    int idx = 0;
    for (int i = 0; i < n; i++) {
        file << xs[i] << " ";
        for (int j = 0; j < stride; j++) {
            idx = i * stride;
            file << yss[idx + j] << " ";
        }
        file << std::endl;
    }

    file.close();
}