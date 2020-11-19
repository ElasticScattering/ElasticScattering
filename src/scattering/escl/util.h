#pragma once

inline double SimpsonWeight(const int i, const int dim) {
    const double main_multiplier = (i % 2 == 0) ? 2.0 : 4.0;
    const bool is_edge = i == 0 || i == (dim - 1);
    return is_edge ? 1.0 : main_multiplier;
}

inline double SimpsonWeight2D(unsigned int i, int j, int dim) {
    return SimpsonWeight(i, dim) * SimpsonWeight(j, dim);
}

inline double GetSigma(double lt, double phi, double tau, double w)
{
    double z = exp(-lt / tau);

    double r = cos(phi) - cos(phi + w * lt) * z;

    r += w * tau * sin(phi + w * lt) * z;
    r -= w * tau * sin(phi);

    return r;
}