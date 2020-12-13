#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <string>
#include "src/scattering/escl/v2.h"


std::vector<v2> GetTestImpurities()
{
    std::filebuf fb;
    if (!fb.open("tests/data/test_impurities.dat", std::ios::in)) {
        printf("Could not load test_impurities.");
        return std::vector<v2>();
    }

    std::vector<v2> impurities;

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
            double x = atof(line.substr(0, pos).c_str());
            double y = atof(line.substr(pos + 1, line.size() - 1).c_str());

            impurities.push_back(v2(x, y));
        }
    }

    fb.close();
    return impurities;
}