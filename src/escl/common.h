#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "src/datastructures/ScatteringParameters.h"
#include "src/datastructures/IntegralResult.h"
#include "constants.h"
#include "weights.h"
#include "device_macros.h"
#include "lifetime.h"

inline double sim_phi_lifetime(const double2 pos, int quadrant, int step, BUFFER_ARGS)
{
    bool clockwise = (sp->is_clockwise == 1);
    bool incoherent = (sp->is_incoherent == 1);
    bool diag_regions = (sp->is_diag_regions == 1);

    const double phi = sp->integrand_start_angle + quadrant * (PI * 0.5) + step * sp->integrand_step_size;

    const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
    return lifetime(min(15.0 * sp->tau, bound_time), pos, phi, clockwise, sp, impurities, cell_index);
}

inline double phi_lifetime(const double2 pos, BUFFER_ARGS)
{
    bool should_compute_sigma = ShouldComputeSigma(sp->mode);
    const bool clockwise = (sp->is_clockwise == 1);
    bool diag_regions = (sp->is_diag_regions == 1);
    bool incoherent = (sp->is_incoherent == 1);
    
    const double w = (clockwise ? -sp->angular_speed : sp->angular_speed);

    double integral = 0;
    for (int j = 0; j < 4; j++)
    {
        const double start = sp->integrand_start_angle + j * (PI * 0.5);
        double total = 0.0;

        for (int i = 0; i < sp->integrand_steps; i++)
        {
            const double phi = start + i * sp->integrand_step_size;

            const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
            double result = lifetime_old(min(15.0 * sp->tau, bound_time), pos, phi, clockwise, sp, impurities, cell_index);

            if (should_compute_sigma) {
                result = GetSigma(result, phi, sp->tau, w, clockwise);

                if (sp->mode != MODE_SIMULATION) {
                    result *= (sp->mode == MODE_SIGMA_XX) ? cos(phi) : sin(phi);
                }
            }

            double w = SimpsonWeight(i, sp->integrand_steps);
            total += result * w;
        }

        integral += total * (sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0));
    }

    return integral;
}

inline double SimpsonWeight(const int i, const int dim) {
    const double main_multiplier = (i % 2 == 0) ? 2.0 : 4.0;
    const bool is_edge = i == 0 || i == (dim - 1);
    return is_edge ? main_multiplier : 1.0;
}

inline double GetWeight2D(unsigned int i, int j, int dim) {
    double w = SimpsonWeight(i, dim) * SimpsonWeight(j, dim);

    return w;
}

inline bool ShouldComputeSigma(int m) {
    return (m == MODE_SIMULATION || m == MODE_SIGMA_XX || m == MODE_SIGMA_XY);
}

inline double GetSigma(double lt, double phi, double tau, double w)
{
    double z = exp(-lt / tau);

    double r = cos(phi) - cos(phi + w * lt) * z;

    r += w * tau * sin(phi + w * lt) * z;
    r -= w * tau * sin(phi);

    return r;
}

#endif // CL_COMMON_H
