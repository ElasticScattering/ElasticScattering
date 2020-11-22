#include "ElasticScattering.h"
#include "escl/lifetime.h"
#include "escl/util.h"

SigmaResult ElasticScatteringCPU::ComputeResult(ScatteringParameters& sp, ImpurityGridIndex& grid)
{
    //CompleteSimulationParameters(sp);

    // GPU kernel works only with even work size.
    const int limit = sp.dim - 1;

    const int values_per_particle = sp.integrand_steps * 4;
    double* matrix = new double[limit * limit * values_per_particle];

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++) {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            for (int q = 0; q < 4; q++) {
                for (int p = 0; p < sp.integrand_steps; p++) {
                    matrix[j * limit + i * values_per_particle + (q * sp.integrand_steps) + p] = lifetime(q, p, pos, &sp, grid.impurities, grid.imp_index);
                }
            }
        }
    }

    // Apply weights for integration.
    const double w = (sp.is_clockwise == 1) ? -sp.angular_speed : sp.angular_speed;

    SigmaResult sr;

    double integral_total = 0;
    for (int j = 0; j < limit; j++) {
        double wy = SimpsonWeight(j, limit);
        for (int i = 0; i < limit; i++) {
            double wx = SimpsonWeight(i, limit);
            for (int q = 0; q < 4; q++) {
                for (int p = 0; p < sp.integrand_steps; p++) {
                    double wp = SimpsonWeight(p, sp.integrand_steps);

                    double lt = matrix[j * limit + i * values_per_particle + (q * sp.integrand_steps) + p];
                    double phi = sp.integrand_start_angle + q * (PI * 0.5) + p * sp.integrand_step_size;

                    double sigma_base = GetSigma(lt, phi, sp.tau, w) * (wy * wx * wp);
                    sr.xx += sigma_base * cos(phi);
                    sr.xy += sigma_base * cos(phi);
                }
            }
        }
    }

    double factor = sp.integrand_angle_area / (4 * (sp.integrand_steps - 1) * (limit * limit));
    factor *= SigmaFactor(sp);
    sr.xx *= factor;
    sr.xy *= factor;

    return sr;
}
