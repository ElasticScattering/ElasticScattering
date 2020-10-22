#ifndef CL_COMMON_H
#define CL_COMMON_H

#include "src/datastructures/ScatteringParameters.h"
#include "constants.h"
#include "weights.h"
#include "device_macros.h"
#include "lifetime.h"

inline double single_lifetime(const double2 pos, const double phi, BUFFER_ARGS) {
    if (sp->angular_speed != 0) {
        const bool clockwise = (sp->is_clockwise == 1);
        const double bound_time = GetBoundTime(phi, sp->alpha, sp->angular_speed, (sp->is_incoherent == 1), (sp->is_diag_regions == 1), clockwise, false);
 
        return lifetimeB(min(15.0 * sp->tau, bound_time), pos, phi, clockwise, sp, impurities, cell_index);
    }
    else {
        return lifetime0(pos, phi, sp, impurities, cell_index);
    }
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

inline double phi_lifetime(const double2 pos, BUFFER_ARGS)
{
    bool should_compute_sigma = ShouldComputeSigma(sp->mode);
    bool clockwise = (sp->is_clockwise == 1);
    double integral = 0;

    for (int j = 0; j < 4; j++) {
        const double start = sp->integrand_start_angle + j * (PI * 0.5);
        double section_total = 0.0;

        for (int i = 0; i < sp->integrand_steps; i++) {
            const double phi = start + i * sp->integrand_step_size;

            double lt = single_lifetime(pos, phi, sp, impurities, cell_index);

            if (should_compute_sigma) {
                lt = GetSigma(lt, phi, sp->tau, sp->angular_speed, clockwise);
                lt *= (sp->mode == MODE_SIGMA_XX) ? cos(phi) : sin(phi);
            }
      
            section_total += lt * GetWeight(i, sp->integrand_steps);
        }

        integral += section_total * sp->integrand_angle_area / ((sp->integrand_steps - 1) * 3.0);
    }

    return integral;
}

#endif // CL_COMMON_H
