#pragma once

#include "ConfigSettings.h"
#include "utils/ErrorMacros.h"

#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <unordered_map>
#include <random>

enum class ProgramMode {
    Test,
    Simulation
};

struct InitParameters
{
    ProgramMode mode;
    std::string config_file;
};

enum class OutputType
{
    Nothing,
    Results,
    All
};

struct SimulationConfiguration
{
    int num_samples;
    int positions_per_row;
    int particles_per_quadrant;
    
    bool print_info;

    OutputType output_type;
    bool use_gpu;

    Settings settings;

    std::vector<double> magnetic_fields;
    std::vector<double> temperatures;

    std::vector<unsigned int> sample_seeds;
    
    std::string config_path;
    std::string base_output_directory;
    std::string output_directory;
    int digits_in_sample_num;

    static SimulationConfiguration ParseFromeFile(std::string file_path)
    {
        std::unordered_map<std::string, std::string> values;

        std::filebuf fb;
        if (!fb.open(file_path, std::ios::in)) {
            std::cout << "Could not load config file '" << file_path << "'";
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

        SimulationConfiguration cfg;
        cfg.config_path = file_path;

        cfg.output_type = OutputType::All;
        if (values.find("output_type") != values.end()) {
            auto log_level = values.at("output_type");
            
            if      (log_level == "nothing") cfg.output_type = OutputType::Nothing;
            else if (log_level == "results") cfg.output_type = OutputType::Results;
        }

        if (cfg.output_type == OutputType::Nothing) 
        {
            cfg.base_output_directory = "";
            cfg.output_directory = "";
        } 
        else 
        {
            cfg.base_output_directory = values.find("output_directory") == values.end() ? "ESLogs" : values.at("output_directory");
            std::string base_dir_name = cfg.base_output_directory + "/Result_";

            unsigned int n = 0;
            while (true)
            {
                std::string dir = base_dir_name + std::to_string(n);
                if (!std::filesystem::exists(dir))
                {
                    cfg.output_directory = dir;
                    break;
                }

                n++;
            }
        }

        bool device_selected = values.find("device") != values.end();
        cfg.use_gpu = device_selected ? values.at("device") == "gpu" : false;

        cfg.print_info = values.find("print_info") != values.end() && atoi(values.at("print_info").c_str()) == 1;

        cfg.num_samples = atoi(values.at("num_samples").c_str());
        CFG_EXIT_CONDITION(cfg.num_samples <= 0, "Samples must be greater than 0.");

        cfg.particles_per_quadrant = atoi(values.at("phi_steps").c_str());
        CFG_EXIT_CONDITION(cfg.particles_per_quadrant % 2 == 0, "Phi steps must be odd.");

        cfg.positions_per_row = atoi(values.at("dimension").c_str());;
        CFG_EXIT_CONDITION(cfg.positions_per_row % 2 == 0, "Dimension must be odd.");

        if (values.find("magnetic_field") != values.end()) {
            cfg.magnetic_fields.push_back(atof(values.at("magnetic_field").c_str()));
        }
        else
        {
            double mf_min       = atof(values.at("magnetic_field_min").c_str());
            double mf_max       = atof(values.at("magnetic_field_max").c_str());
            double mf_n         = atoi(values.at("magnetic_field_n").c_str());

            CFG_EXIT_CONDITION(mf_n < 2, "Tried to declare magnetic field range with n<2. Use 'magnetic_field x' when using only one value.");
            cfg.magnetic_fields.resize(mf_n);

            double mf_step_size = (mf_max - mf_min) / (double)(mf_n - 1);
            
            if (cfg.print_info) std::cout << "Using magnetic fields: [";
            for (int i = 0; i < mf_n; i++) {
                cfg.magnetic_fields[i] = mf_min + i * mf_step_size;
                if (cfg.print_info) 
                {
                    std::cout << cfg.magnetic_fields[i];
                    if (i != (mf_n - 1)) {
                        std::cout << ", ";
                    }
                }
            }

            if (cfg.print_info) std::cout << "]." << std::endl;
        }

        {
            std::stringstream string_stream(values.at("temperatures"));

            while (string_stream.good())
            {
                std::string a;
                std::getline(string_stream, a, ' ');
                double temp = atof(a.c_str());
                if (temp > 0)
                    cfg.temperatures.push_back(temp);
            }
        }

        cfg.settings.particle_speed         = atof(values.at("particle_speed").c_str());
        cfg.settings.tau                    = atof(values.at("tau").c_str());
        cfg.settings.alpha                  = atof(values.at("alpha").c_str());
        cfg.settings.is_clockwise           = atoi(values.at("clockwise").c_str()) > 0;

        std::random_device random_device;
        cfg.sample_seeds.resize(cfg.num_samples);
        unsigned int seed = values.find("start_seed") != values.end() ? atoi(values.at("start_seed").c_str()) : random_device();
        if (cfg.print_info) std::cout << "RNG seed: " << seed << std::endl;

        for (int i = 0; i < cfg.num_samples; i++)
            cfg.sample_seeds[i] = (seed+100) * (i+1);

        cfg.settings.region_size            = atof(values.at("region_size").c_str());
        cfg.settings.region_extends         = atof(values.at("region_extends").c_str());
        cfg.settings.impurity_density       = atof(values.at("impurity_density").c_str());
        cfg.settings.impurity_radius        = atof(values.at("impurity_radius").c_str());
        cfg.settings.target_cell_population = values.find("target_cell_population") != values.end() ? atoi(values.at("target_cell_population").c_str()) : 15;

        int x = cfg.num_samples;
        cfg.digits_in_sample_num = 0;
        while (x > 0)
        {
            x = floor(x / 10);
            cfg.digits_in_sample_num++;
        }
        
        std::cout << std::endl;

        return cfg;
    }
};
