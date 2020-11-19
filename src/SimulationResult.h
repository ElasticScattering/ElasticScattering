#pragma once

#include <vector>

struct SigmaResult {
    double xx;
    double xy;

    SigmaResult() {
        xx = 0;
        xy = 0;
    }
};

struct TextureResults
{
    std::vector<double> xx, xy, xx_inc, xy_inc;

    TextureResults(int size) {
        xx.resize(size); 
        xy.resize(size); 
        xx_inc.resize(size);
        xy_inc.resize(size);
    }
};
