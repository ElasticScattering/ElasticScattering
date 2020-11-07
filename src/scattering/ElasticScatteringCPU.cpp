#include "ElasticScattering.h"
#include "escl/lifetime.h"
#include "escl/util.h"

SigmaResult ElasticScatteringCPU::ComputeResult(ScatteringParameters& p_sp)
{
    PrepareCompute(p_sp);

    // GPU kernel works only with even work size.
    const int limit = sp.dim - 1;
    
    const int values_per_particle = p_sp.integrand_steps * 4;
    double* matrix = new double[limit * limit * values_per_particle];

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++) {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            for (int q = 0; q < 4; q++) {
                for (int p = 0; p < p_sp.integrand_steps; p++) {
                    matrix[j * limit + i * values_per_particle + (q*p_sp.integrand_steps) + p] = sim_phi_lifetime(pos, q, p, &p_sp, grid.impurities, grid.imp_index);
                }
            }
        }
    }

    // Apply weights for integration.
    const double w = (p_sp.is_clockwise == 1) ? -p_sp.angular_speed : p_sp.angular_speed;

    SigmaResult sr;

    double integral_total = 0;
    for (int j = 0; j < limit; j++) {
        double wy = SimpsonWeight(j, limit);
        for (int i = 0; i < limit; i++) {
            double wx = SimpsonWeight(i, limit);
            for (int q = 0; q < 4; q++) {
                for (int p = 0; p < p_sp.integrand_steps; p++) {
                    double wp = SimpsonWeight(p, p_sp.integrand_steps);

                    double lt = matrix[j * limit + i * values_per_particle + (q * p_sp.integrand_steps) + p];
                    double phi = p_sp.integrand_start_angle + q * (PI * 0.5) + p * p_sp.integrand_step_size;

                    double sigma_base = GetSigma(lt, phi, p_sp.tau, w) * (wy * wx * wp);
                    sr.xx += sigma_base * cos(phi);
                    sr.xy += sigma_base * cos(phi);
                }
            }
        }
    }

    double factor = p_sp.integrand_angle_area / (4 * (p_sp.integrand_steps - 1) * (limit * limit));
    factor *= SigmaFactor();
    sr.xx *= factor;
    sr.xy *= factor;

    return sr;
}

bool ElasticScatteringCPU::PrepareCompute(ScatteringParameters &p_sp) {
    CompleteSimulationParameters(p_sp);
    
    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed  = (sp.dim != p_sp.dim);

    if (!first_run && !impurities_changed && !work_size_changed) return false;

    sp = p_sp;
        
    if (first_run || impurities_changed)
        grid.Generate(sp);

    first_run = false;

    return true;
}
