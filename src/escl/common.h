#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "src/datastructures/ScatteringParameters.h"
#include "src/datastructures/IntegralResult.h"
#include "constants.h"
#include "weights.h"
#include "device_macros.h"
#include "lifetime.h"

inline ScatterResult sim_phi_lifetime(const double2 pos, BUFFER_ARGS)
{
    bool should_compute_sigma = ShouldComputeSigma(sp->mode);
    bool clockwise = (sp->is_clockwise == 1);
    bool incoherent = (sp->is_incoherent == 1);
    bool diag_regions = (sp->is_diag_regions == 1);
    
    ScatterResult result;
    result.xx = 0;
    result.xy = 0;

    for (int j = 0; j < 4; j++) {
        const double start = sp->integrand_start_angle + j * (PI * 0.5);
        double2 section_total = { 0.0, 0.0 };

        for (int i = 0; i < sp->integrand_steps; i++) {
            const double phi = start + i * sp->integrand_step_size;

            const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, incoherent, diag_regions, clockwise, false);
            const double lt = lifetime_old(min(15.0 * sp->tau, bound_time), pos, phi, clockwise, sp, impurities, cell_index);

            double result = GetSigma(lt, phi, sp->tau, sp->angular_speed, clockwise);
            double xx = result * cos(phi);
            double xy = result * sin(phi);
      
            double w = GetWeight(i, sp->integrand_steps);
            section_total = section_total + (double2)(xx, xy) * w;
        }

        section_total = section_total * (sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0));

        result.xx += section_total.x;
        result.xy += section_total.y;
    }

    return result;
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

            double w = GetWeight(i, sp->integrand_steps);
            total += result * w;
        }

        integral += total * (sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0));
    }

    return integral;
}

inline bool ShouldComputeSigma(int m) {
    return (m == MODE_SIMULATION || m == MODE_SIGMA_XX || m == MODE_SIGMA_XY);
}

inline double GetSigma(double lt, double phi, double tau, double w, bool clockwise)
{
    w = clockwise ? -w : w;

    double z = exp(-lt / tau);

    double r = cos(phi) - cos(phi + w * lt) * z;

    r += w * tau * sin(phi + w * lt) * z;
    r -= w * tau * sin(phi);

    return r;
}

#endif // CL_COMMON_H
