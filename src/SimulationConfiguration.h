#pragma once

#include "Settings.h"

#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <unordered_map>


enum class ProgramMode {
    Test,
    Simulation
};

typedef struct
{
    ProgramMode mode;
    bool use_gpu;
    bool dont_show_info;
    bool write_images;
    std::string config_file;
} InitParameters;

typedef struct Range
{
    double min, max, step_size;
    int n;
} Range;

typedef struct SimulationConfiguration
{
    int num_samples;
    std::vector<double> magnetic_fields;
    std::vector<double> temperatures;

    int particles_per_row;
    int quadrant_integral_steps;

    Settings settings;

    bool log_intermediates;
    std::string base_output_directory;
    std::string output_directory;
    int digits_in_sample_num;

    static SimulationConfiguration ParseFromeFile(std::string file_path)
    {
        std::unordered_map<std::string, std::string> values;

        std::filebuf fb;
        if (!fb.open(file_path, std::ios::in)) {
            std::cout << "Could not load config file.";
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

        cfg.log_intermediates = atoi(values.at("log_intermediates").c_str()) == 1;
        cfg.num_samples = atoi(values.at("num_samples").c_str());
        cfg.quadrant_integral_steps = atoi(values.at("integrand_steps").c_str());
        cfg.particles_per_row = atoi(values.at("dimension").c_str());;

        {
            Range magnetic_field;
            magnetic_field.min       = atof(values.at("magnetic_field_min").c_str());
            magnetic_field.max       = atof(values.at("magnetic_field_max").c_str());
            magnetic_field.n         = atoi(values.at("magnetic_field_n").c_str());
            magnetic_field.step_size = (magnetic_field.max - magnetic_field.min) / (double)(magnetic_field.n - 1);

            cfg.magnetic_fields.resize(magnetic_field.n);
            for (int i = 0; i < magnetic_field.n; i++) {
                cfg.magnetic_fields[i] = magnetic_field.min + i * magnetic_field.step_size;
            }
        }

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

        cfg.settings.particle_speed                  = atof(values.at("particle_speed").c_str());
        cfg.settings.tau                             = atof(values.at("tau").c_str());
        cfg.settings.alpha                           = atof(values.at("alpha").c_str());
        cfg.settings.is_clockwise                    = atoi(values.at("clockwise").c_str()) > 0;

        cfg.settings.region_size                     = atof(values.at("region_size").c_str());
        cfg.settings.region_extends                  = atof(values.at("region_extends").c_str());
        cfg.settings.impurity_density                = atof(values.at("impurity_density").c_str());
        cfg.settings.impurity_radius                 = atof(values.at("impurity_radius").c_str());
        cfg.settings.max_expected_impurities_in_cell = atoi(values.at("max_expected_impurities_in_cell").c_str());

        int x = cfg.num_samples;
        cfg.digits_in_sample_num = 0;
        while (x > 0)
        {
            x = floor(x / 10);
            cfg.digits_in_sample_num++;
        }

        return cfg;
    }
} SimulationConfiguration;
