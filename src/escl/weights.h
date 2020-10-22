#pragma once

inline bool IsEdge(int i, int dim) {
    return (i == 0) || (i == (dim - 1));
}

inline double GetWeightValue(int i) {
    return ((i % 2) == 0) ? 2.0 : 4.0;
}

inline double GetWeight(int i, int dim) {
    bool is_edge = IsEdge(i, dim);
    double main_multiplier = GetWeightValue(i);
    double w = is_edge ? main_multiplier : 1.0;

    return w;
}

inline double GetWeight2D(unsigned int i, int j, int dim) {
    double w = GetWeight(i, dim) * GetWeight(j, dim);
    
    return w;
}