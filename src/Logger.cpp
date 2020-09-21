#include "src/Logger.h"
#include "src/escl/constants.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <iomanip> 
#include <ctime>

namespace fs = std::filesystem;

void Logger::LogResult(const SimulationParameters& sp, const SimulationResult& sr) const
{
    auto date_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::filesystem::create_directory("ES Logs");
    
    std::ofstream file;

    for (unsigned int n = 0; ; ++n) { //todo leidende nullen 0008
        std::string fname = std::string("Es Logs/results_") + std::to_string(n) + std::string(".dat");

        std::ifstream ifile;
        ifile.open(fname.c_str());

        if (!ifile.is_open()) {
            file.open(fname.c_str());
            break;
        }

        ifile.close();
    }

    file << std::scientific << std::setprecision(10);
    file << "# Elastic Scattering Results" << std::endl;
    file << "# Completed on: " << std::ctime(&date_time) << "." << std::endl;
    file << "# Elapsed time: " << sr.time_elapsed << " seconds." << std::endl;
    file << "# Each row is the average of "<< sr.iterations_per_run << " iterations with different impurity seeds): " << std::endl;
    file << "# Simulation parameters:" << std::endl;
    file << "#\t" << "Integrand steps:  " << sp.integrand_steps << std::endl;
    file << "#\t" << "Dimension:        " << sp.dim << std::endl;
    file << "#\t" << "Diag. regions:    " << ((sp.is_diag_regions == 1) ? "True" : "False") << std::endl;
    file << "#\t" << "Clockwise:        " << ((sp.is_clockwise == 1)    ? "True" : "False") << std::endl;
    file << "#\t" << "Incoherent:       " << ((sp.is_incoherent == 1)   ? "True" : "False") << std::endl; 
    file << "#" << std::endl;
    file << "#\t" << "Temperature:      " << sp.temperature << std::endl;
    file << "#\t" << "Tau:              " << sp.tau << std::endl;
    file << "#\t" << "Magnetic field:   " << sp.magnetic_field << std::endl;
    file << "#\t" << "Phi:              " << sp.phi << std::endl;
    file << "#\t" << "Alpha:            " << sp.alpha << std::endl;
    file << "#\t" << "Particle speed:   " << sp.particle_speed << std::endl;
    file << "#\t" << "Angular speed:    " << sp.angular_speed << std::endl;
    file << "#\n# Impurities:" << std::endl;
    file << "#\t" << "Region size:      " << sp.region_size << std::endl;
    file << "#\t" << "Region extends:   " << sp.region_extends << std::endl; 
    file << "#\t" << "Density:          " << sp.impurity_density << std::endl;
    file << "#\t" << "Count:            " << sp.impurity_count << std::endl;
    file << "#\t" << "Radius:           " << sp.impurity_radius << std::endl;

    file << "#\n#Constants:" << std::endl;
    file << "#\t" << "Mass0: " << M0 << std::endl;
    file << "#\t" << "E:     " << E << std::endl;
    file << "#\t" << "HBAR:  " << HBAR << std::endl;
    file << "#\t" << "C:     " << C1 << std::endl;
    file << "#\t" << "KB:    " << KB << std::endl;

    file << "#\n#Results:\n" << std::endl;

    file << "magnetic_field temperature sigma_xx_inc sigma_xx_coh sigma_xy_inc sigma_xy_coh delta_xx" << std::endl;

    int n = sr.xs.size();

    int idx = 0;
    for (int i = 0; i < n; i++) {
        file << sr.xs[i] << " " << sr.xs_temperature[i];
        file << " " << sr.results_xxi[i] << " " << sr.results_xx[i];
        file << " " << sr.results_xyi[i] << " " << sr.results_xy[i];
        file << " " << sr.delta_xxi[i] << std::endl;
    }

    file.flush();
    file.close();
}