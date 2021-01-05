#include "windows.h"
#include "Simulation.h"
#include "src/utils/OpenCLUtils.h"

#include <iostream> // @Todo, remove.
#include <iomanip>

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_command_queue queue;
} OCLResources;

typedef struct
{
    cl_program program_lifetimes;
    cl_kernel lifetimes_kernel;

    cl_mem raw_lifetimes;

    cl_mem particle_settings;
    cl_mem impurity_settings;
    cl_mem simulation_settings;

    cl_mem impurities;
    cl_mem cell_indices;
    cl_mem metrics;
} OCLSimResources;

typedef struct
{
    cl_program program_integration;

    cl_kernel apply_max_lifetime;
    cl_kernel apply_sigma_comp_kernel;
    cl_kernel apply_simpson_weights;
    cl_kernel integrate_position_kernel;
    cl_kernel sum_kernel;

    cl_mem lifetimes;

    cl_mem lifetimes_sigma_xx;
    cl_mem lifetimes_sigma_xy;

    cl_mem lifetimes_positions;
    cl_mem sigma_xx_positions;
    cl_mem sigma_xy_positions;

    cl_mem sigma_xx;
    cl_mem sigma_xy;

    cl_mem incomplete_sum; // Nog steeds nodig?
} OCLIntegrationResources;

OCLResources ocl;
OCLSimResources ocl_scatter;
OCLIntegrationResources ocl_integration;


std::vector<Sigma> SimulationCL::ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
    Metrics metrics;
    QueryPerformanceCounter(&pc.lifetimeBegin);
    ComputeLifetimes(magnetic_field, grid, metrics);
    QueryPerformanceCounter(&pc.lifetimeEnd);
    metrics.time_elapsed_lifetimes = GetElapsedTime(pc.lifetimeBegin, pc.lifetimeEnd);

    metrics.real_particles = ss.total_particles - metrics.particle_metrics.particles_inside_impurity;

    QueryPerformanceCounter(&pc.temperaturesBegin);
    std::vector<Sigma> results(temperatures.size());
    for (int i = 0; i < temperatures.size(); i++)
        results[i] = DeriveTemperature(temperatures[i]);
    QueryPerformanceCounter(&pc.temperaturesEnd);
    metrics.time_elapsed_temperatures = GetElapsedTime(pc.temperaturesBegin, pc.temperaturesEnd);

    sample_metrics.iteration_metrics.push_back(metrics);

    return results;
}

std::vector<IterationResult> SimulationCL::ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
    Metrics metrics;
    QueryPerformanceCounter(&pc.lifetimeBegin);
    ComputeLifetimes(magnetic_field, grid, metrics);
    QueryPerformanceCounter(&pc.lifetimeEnd);
    metrics.time_elapsed_lifetimes = GetElapsedTime(pc.lifetimeBegin, pc.lifetimeEnd);

    metrics.real_particles = ss.total_particles - metrics.particle_metrics.particles_inside_impurity;

    QueryPerformanceCounter(&pc.temperaturesBegin);
    std::vector<IterationResult> results(temperatures.size());
    for (int i = 0; i < temperatures.size(); i++)
        results[i] = DeriveTemperatureWithImages(temperatures[i]);
    QueryPerformanceCounter(&pc.temperaturesEnd);
    metrics.time_elapsed_temperatures = GetElapsedTime(pc.temperaturesBegin, pc.temperaturesEnd);

    sample_metrics.iteration_metrics.push_back(metrics);

    return results;
}

void SimulationCL::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
    ps.angular_speed = E * magnetic_field / M;
    ss.signed_angular_speed = ps.is_clockwise ? -ps.angular_speed : ps.angular_speed;

    if (grid.GetSeed() != last_grid_seed) {
        UploadImpurities(grid);
        last_grid_seed = grid.GetSeed();
    }

    cl_int clStatus;
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.particle_settings,   CL_TRUE, 0, sizeof(ParticleSettings),   (void*)&ps, 0, nullptr, nullptr);
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.impurity_settings,   CL_TRUE, 0, sizeof(ImpuritySettings),   (void*)&is, 0, nullptr, nullptr);
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.simulation_settings, CL_TRUE, 0, sizeof(SimulationSettings), (void*)&ss, 0, nullptr, nullptr);

    ParticleMetrics pm;
    ocl_scatter.metrics = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(SimulationSettings), nullptr, &clStatus);
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.metrics, CL_TRUE, 0, sizeof(ParticleMetrics), (void*)&pm, 0, nullptr, nullptr);

    // Ipv WI te maken die elk een positie behandelt, laten we elk WI één particle doen.
    // @Optimize Verifiëer dat dit ook daadwerkelijk achtereenvolgende phi's berekend in één WG.
    {
        work_size.positions_per_row      = (size_t)ss.positions_per_row + 1;
        work_size.particles_per_position = (size_t)ss.particles_per_position;
        work_size.total_particles        = work_size.positions_per_row * work_size.positions_per_row * work_size.particles_per_position;

        work_size.particles_global[0] = work_size.positions_per_row;
        work_size.particles_global[1] = work_size.positions_per_row;
        work_size.particles_global[2] = work_size.particles_per_position;

        const size_t values_in_work_group = min(work_size.particles_per_position, 256); // @Todo, dit kan anders zijn voor andere devices.

        size_t other = 256 / values_in_work_group;
        work_size.particles_local[0] = max(1, (other % 2 == 0) ? other : other - 1);
        work_size.particles_local[1] = 1;
        work_size.particles_local[2] = values_in_work_group;

        std::cout << "Self determined work size:" << std::endl;
        std::cout << "Global: " << work_size.particles_global[0] << " " << work_size.particles_global[1] << " " << work_size.particles_global[2] << std::endl;
        std::cout << "Local:  " << work_size.particles_local[0]  << " " << work_size.particles_local[1]  << " " << work_size.particles_local[2]  << std::endl;
    }

    ocl_scatter.raw_lifetimes = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create metrics buffer.");

    {
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 0, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 1, sizeof(cl_mem), (void*)&ocl_scatter.particle_settings);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 2, sizeof(cl_mem), (void*)&ocl_scatter.impurity_settings);

        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 3, sizeof(cl_mem), (void*)&ocl_scatter.impurities);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 4, sizeof(cl_mem), (void*)&ocl_scatter.cell_indices);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 5, sizeof(cl_mem), (void*)&ocl_scatter.raw_lifetimes);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 6, sizeof(cl_mem), (void*)&ocl_scatter.metrics);
    }

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_scatter.lifetimes_kernel, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start lifetimes kernel.");
    clFinish(ocl.queue);

    //std::vector<ParticleMetrics> metrics_holder(1);
    clEnqueueReadBuffer(ocl.queue, ocl_scatter.metrics, CL_TRUE, 0, sizeof(ParticleMetrics), &metrics.particle_metrics, 0, nullptr, nullptr);
    //metrics.particle_metrics = metrics_holder[0];

    //@Todo: remove?
    std::vector<double> lifetimes(work_size.total_particles);
    clStatus = clEnqueueReadBuffer(ocl.queue, ocl_scatter.raw_lifetimes, CL_TRUE, 0, sizeof(double) * work_size.total_particles, lifetimes.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");

    std::ofstream file;
    file.open("verification/gpu_lt.txt");
    file << std::setprecision(12);
    for (int j = 0; j < ss.positions_per_row; j++) {
        file << std::endl << "R" << j << std::endl;
        for (int i = 0; i < ss.positions_per_row; i++) {
            file << std::endl;
            for (int q = 0; q < 4; q++) {
                for (int p = 0; p < ss.particles_per_quadrant; p++) {
                    int idx = j * work_size.positions_per_row * ss.particles_per_position + i * ss.particles_per_position + q * ss.particles_per_quadrant + p;
                    
                    file << lifetimes[idx] << std::endl;
                }
            }
        }
    }

    metrics.avg_particle_lifetime = AverageLifetime(lifetimes);
}

Sigma SimulationCL::DeriveTemperature(const double temperature) const
{
    cl_int clStatus;
    double tau = GetTau(temperature);
    
    ocl_integration.lifetimes = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create integration lifetimes buffer.");

    ocl_integration.lifetimes_sigma_xx = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sigma_lifetimes_xx buffer.");

    ocl_integration.lifetimes_sigma_xy = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sigma_lifetimes_xy buffer.");

    // Apply max lifetime.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 0, sizeof(cl_mem), (void*)&ocl_scatter.raw_lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 1, sizeof(double), (void*)&ocl_scatter.simulation_settings);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        double default_max_lifetime = GetDefaultMaxLifetime(tau);
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 2, sizeof(double), (void*)&default_max_lifetime);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 3, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_max_lifetime, 3, nullptr, work_size.particles_global, work_size.particles_local, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start apply_max_lifetime execution.");

#ifdef _DEBUG
        clFinish(ocl.queue);
        std::vector<double> truncated_lifetimes(work_size.total_particles);
        clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes, CL_TRUE, 0, sizeof(double) * work_size.total_particles, truncated_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
#endif
    }

    // Apply sigma components.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 2, sizeof(cl_mem), (void*)&ocl_scatter.particle_settings);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 3, sizeof(double), (void*)&tau);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        int mode = 0;

        {
            mode = MODE_SIGMA_XX;
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 4, sizeof(int), (void*)&mode);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 5, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xx);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_sigma_comp_kernel, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_sigma_comp_kernel xx execution.");
        }

        {
            mode = MODE_SIGMA_XY;
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 4, sizeof(int), (void*)&mode);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 5, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xy);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_sigma_comp_kernel, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_sigma_comp_kernel xy execution.");
        }

#ifdef _DEBUG
        clFinish(ocl.queue);
        std::vector<double> sigma_xx_lifetimes(work_size.total_particles);
        clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xx, CL_TRUE, 0, sizeof(double) * work_size.total_particles, sigma_xx_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");

        std::vector<double> sigma_xy_lifetimes(work_size.total_particles);
        clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xy, CL_TRUE, 0, sizeof(double) * work_size.total_particles, sigma_xy_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
#endif
    }

    // Apply simpson weights.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 2, sizeof(int), (void*)&ss.particles_per_quadrant);

        {
            clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xx);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_simpson_weights, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_simpson_weights execution.");
        }

        {
            clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xy);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_simpson_weights, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_simpson_weights execution.");
        }
    }

    // @Temporary, gebruik sum kernel.
    clFinish(ocl.queue);
    std::vector<double> sigma_xx_lifetimes_weighted(work_size.total_particles);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xx, CL_TRUE, 0, sizeof(double)* work_size.total_particles, sigma_xx_lifetimes_weighted.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");

    std::vector<double> sigma_xy_lifetimes_weighted(work_size.total_particles);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xy, CL_TRUE, 0, sizeof(double)* work_size.total_particles, sigma_xy_lifetimes_weighted.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");

    Sigma sigma;
    for (int i = 0; i < sigma_xx_lifetimes_weighted.size(); i++)
        sigma.xx += sigma_xx_lifetimes_weighted[i];

    for (int i = 0; i < sigma_xy_lifetimes_weighted.size(); i++)
        sigma.xy += sigma_xy_lifetimes_weighted[i];

    double factor = GetSigmaIntegrandFactor(tau);
    sigma.xx *= factor;
    sigma.xy *= factor;

    return sigma;

    //Sigma sigma;
    // Sum sigma xx/xy buffers.
    // @Optimize, voer de som nog vaker uit tot een aantal elementen.
    /*
    {
        const size_t half_size = work_size.total_particles / 2;
        const size_t incomplete_sum_size = half_size / work_size.particles_per_position; //Items in een WI....
        ocl_integration.incomplete_sum = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * incomplete_sum_size, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create metrics buffer.");

        clStatus = clSetKernelArg(ocl_integration.sum_kernel, 1, sizeof(double) * work_size.particles_per_position, nullptr);
        clStatus = clSetKernelArg(ocl_integration.sum_kernel, 2, sizeof(cl_mem), (void*)&ocl_integration.incomplete_sum);
        
        {
            clStatus = clSetKernelArg(ocl_integration.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xx);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.sum_kernel, 1, nullptr, &half_size, &work_size.particles_per_position, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start sum_kernel execution.");

            std::vector<double> incomplete_sum(incomplete_sum_size);
            clEnqueueReadBuffer(ocl.queue, ocl_integration.incomplete_sum, CL_TRUE, 0, sizeof(double) * incomplete_sum_size, incomplete_sum.data(), 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Failed to read back incomplete sum XX.");

            for (int i = 0; i < incomplete_sum.size(); i++)
                sigma.xx += incomplete_sum[i];
        }

        //clEnqueueFillBuffer(ocl.queue, ocl_scatter.incomplete_sum, 0, 0, 0, sizeof(double) * work_size.total_particles, 0, 0, 0);

        {
            clStatus = clSetKernelArg(ocl_integration.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xy);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.sum_kernel, 1, nullptr, &half_size, &work_size.particles_per_position, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start sum_kernel execution.");

            std::vector<double> incomplete_sum(incomplete_sum_size);
            clEnqueueReadBuffer(ocl.queue, ocl_integration.incomplete_sum, CL_TRUE, 0, sizeof(double) * incomplete_sum_size, incomplete_sum.data(), 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Failed to read back incomplete sum XY.");

            for (int i = 0; i < incomplete_sum.size(); i++)
                sigma.xy += incomplete_sum[i];
        }
    }

    double factor = GetSigmaIntegrandFactor(tau);
    sigma.xx *= factor;
    sigma.xy *= factor;

    return sigma;
    */

}

IterationResult SimulationCL::DeriveTemperatureWithImages(const double temperature) const
{
    cl_int clStatus;
    double tau = GetTau(temperature);

    ocl_integration.lifetimes = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create integration lifetimes buffer.");

    ocl_integration.lifetimes_sigma_xx = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sigma_lifetimes_xx buffer.");

    ocl_integration.lifetimes_sigma_xy = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * work_size.total_particles, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sigma_lifetimes_xy buffer.");

    int total_positions = work_size.positions_per_row * work_size.positions_per_row;

    ocl_integration.lifetimes_positions = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * total_positions, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes_particle buffer.");

    ocl_integration.sigma_xx_positions = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * total_positions, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sigma_lifetimes_xy buffer.");

    ocl_integration.sigma_xy_positions = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * total_positions, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sigma_lifetimes_xy buffer.");

    // Apply max lifetime.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 0, sizeof(cl_mem), (void*)&ocl_scatter.raw_lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 1, sizeof(double), (void*)&ocl_scatter.simulation_settings);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        double default_max_lifetime = GetDefaultMaxLifetime(tau);
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 2, sizeof(double), (void*)&default_max_lifetime);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 3, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_max_lifetime, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start apply_max_lifetime execution.");

#ifdef _DEBUG
        clFinish(ocl.queue);
        std::vector<double> truncated_lifetimes(work_size.total_particles);
        clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes, CL_TRUE, 0, sizeof(double) * work_size.total_particles, truncated_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
#endif
    }

    // Apply sigma components.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 2, sizeof(cl_mem), (void*)&ocl_scatter.particle_settings);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 3, sizeof(double), (void*)&tau);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

        int mode = 0;

        {
            mode = MODE_SIGMA_XX;
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 4, sizeof(int), (void*)&mode);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 5, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xx);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_sigma_comp_kernel, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_sigma_comp_kernel xx execution.");
        }

        {
            mode = MODE_SIGMA_XY;
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 4, sizeof(int), (void*)&mode);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 5, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xy);
            CL_FAIL_CONDITION(clStatus, "Couldn't set argument to kernel.");

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_sigma_comp_kernel, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_sigma_comp_kernel xy execution.");
        }

#ifdef _DEBUG
        clFinish(ocl.queue);
        std::vector<double> sigma_xx_lifetimes(work_size.total_particles);
        clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xx, CL_TRUE, 0, sizeof(double) * work_size.total_particles, sigma_xx_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");

        std::vector<double> sigma_xy_lifetimes(work_size.total_particles);
        clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xy, CL_TRUE, 0, sizeof(double) * work_size.total_particles, sigma_xy_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
#endif
    }

    // Apply simpson weights.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 2, sizeof(int), (void*)&ss.particles_per_quadrant);

        {
            clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xx);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_simpson_weights, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_simpson_weights execution.");
        }

        {
            clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xy);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_simpson_weights, 3, nullptr, work_size.particles_global, NULL, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_simpson_weights execution.");
        }
    }

    {
        // Integrate to position for images.
        clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);

        size_t p = min(work_size.positions_per_row, 256);
        size_t positions_global[2] = { work_size.positions_per_row, work_size.positions_per_row };
        size_t positions_local[2] = { p, max(1, 256/p) };
        
        
        {
            clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);
            clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 2, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_positions);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.integrate_position_kernel, 2, nullptr, positions_global, positions_local, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start integrate_position_kernel execution.");
        } 
        
        {
            clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xx);
            clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 2, sizeof(cl_mem), (void*)&ocl_integration.sigma_xx_positions);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.integrate_position_kernel, 2, nullptr, positions_global, positions_local, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start integrate_position_kernel execution.");
        }

        {
            clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes_sigma_xy);
            clStatus = clSetKernelArg(ocl_integration.integrate_position_kernel, 2, sizeof(cl_mem), (void*)&ocl_integration.sigma_xy_positions);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.integrate_position_kernel, 2, nullptr, positions_global, positions_local, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start integrate_position_kernel execution.");
        }
    }

    clFinish(ocl.queue);

    IterationResult ir;
    std::vector<double> lifetimes_positions(total_positions);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_positions, CL_TRUE, 0, sizeof(double) * total_positions, lifetimes_positions.data(), 0, nullptr, nullptr);
    
    std::vector<double> lifetime_image(ss.total_positions);
    for (int j = 0; j < ss.positions_per_row; j++)
        for (int i = 0; i < ss.positions_per_row; i++) {
            lifetime_image[j * ss.positions_per_row + i] = lifetimes_positions[j * work_size.positions_per_row + i];
        }
    
    ir.particle_lifetimes = lifetime_image;

    std::vector<double> sigma_xx_positions(total_positions);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.sigma_xx_positions, CL_TRUE, 0, sizeof(double) * total_positions, sigma_xx_positions.data(), 0, nullptr, nullptr);
    std::vector<double> xx_image(ss.total_positions);
    for (int j = 0; j < ss.positions_per_row; j++)
        for (int i = 0; i < ss.positions_per_row; i++)
            xx_image[j * ss.positions_per_row + i] = sigma_xx_positions[j * work_size.positions_per_row + i];

    ir.sigmas.xx_buffer = xx_image;

    std::vector<double> sigma_xy_positions(total_positions);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.sigma_xy_positions, CL_TRUE, 0, sizeof(double) * total_positions, sigma_xy_positions.data(), 0, nullptr, nullptr);

    std::vector<double> xy_image(ss.total_positions);
    for (int j = 0; j < ss.positions_per_row; j++)
        for (int i = 0; i < ss.positions_per_row; i++)
            xy_image[j * ss.positions_per_row + i] = sigma_xy_positions[j * work_size.positions_per_row + i];

    ir.sigmas.xy_buffer = xy_image;

    // @Temporary, gebruik sum kernel.
    std::vector<double> sigma_xx_lifetimes_weighted(work_size.total_particles);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xx, CL_TRUE, 0, sizeof(double) * work_size.total_particles, sigma_xx_lifetimes_weighted.data(), 0, nullptr, nullptr);

    std::vector<double> sigma_xy_lifetimes_weighted(work_size.total_particles);
    clEnqueueReadBuffer(ocl.queue, ocl_integration.lifetimes_sigma_xy, CL_TRUE, 0, sizeof(double) * work_size.total_particles, sigma_xy_lifetimes_weighted.data(), 0, nullptr, nullptr);

    for (int i = 0; i < sigma_xx_lifetimes_weighted.size(); i++)
        ir.result.xx += sigma_xx_lifetimes_weighted[i];

    for (int i = 0; i < sigma_xy_lifetimes_weighted.size(); i++)
        ir.result.xy += sigma_xy_lifetimes_weighted[i];

    double factor = GetSigmaIntegrandFactor(tau);
    ir.result.xx *= factor;
    ir.result.xy *= factor;
    
    return ir;
}

void SimulationCL::UploadImpurities(const Grid& grid)
{
    cl_int clStatus;

    clReleaseMemObject(ocl_scatter.impurities);
    clReleaseMemObject(ocl_scatter.cell_indices);

    ocl_scatter.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(v2) * grid.GetTotalImpurityCount(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create impurities buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.impurities, CL_TRUE, 0, sizeof(v2) * grid.GetTotalImpurityCount(), grid.GetImpurities().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't upload impurities.");

    ocl_scatter.cell_indices = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(int) * grid.GetIndex().size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create cell index buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.cell_indices, CL_TRUE, 0, sizeof(int) * grid.GetIndex().size(), grid.GetIndex().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't upload cell indices.");
}

SimulationCL::SimulationCL(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant)
{
    InitializeOpenCL(true, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (false)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    LARGE_INTEGER compileBegin, compileEnd;
    QueryPerformanceCounter(&compileBegin);
    std::cout << "\nBuilding main program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "src/sim/cl/scatter.cl", &ocl_scatter.program_lifetimes);
    QueryPerformanceCounter(&compileEnd);
    double elapsed = GetElapsedTime(compileBegin, compileEnd);
    std::cout << "Compilation succeeded! (" << elapsed << "s)" << std::endl;

    QueryPerformanceCounter(&compileBegin);
    std::cout << "\nBuilding integration program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "src/sim/cl/integration.cl", &ocl_integration.program_integration);
    QueryPerformanceCounter(&compileEnd);
    elapsed = GetElapsedTime(compileBegin, compileEnd);
    std::cout << "Compilation succeeded! (" << elapsed << "s)" << std::endl;

    /*
    //auto clStatus = clBuildProgram(program, 1, &devices[0], "-cl-mad-enable", NULL, NULL); //build the program
    //auto clStatus = clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_NUM_DEVICES, sizeof(size_t), &nb_devices, &nbread);// Return 1 devices

    size_t* np = new size_t[1]; //Create size array
    auto clStatus = clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), np, NULL);//Load in np the size of my binary
    char** bn = new char* [1]; //Create the binary array
    for (int i = 0; i < 1; i++)  bn[i] = new char[np[i]]; // I know... it's bad... but if i use new char[np], i have a segfault... :/
    clStatus = clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARIES, sizeof(unsigned char*), bn, NULL); //Load the binary itself

    printf("%s\n", bn[0]); //Print the first binary. But here, I have some curious characters
    FILE* fp;
    fopen_s(&fp, "binar.bin", "wb");
    fwrite(bn[0], sizeof(char), np[0], fp); // Save the binary, but my file stay empty
    */

    /*
    size_t number_of_binaries;
    clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &number_of_binaries, NULL);
    size_t* binaries = new size_t[number_of_binaries];
    
    auto binary = new unsigned char *[number_of_binaries];
    clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARIES, number_of_binaries, &binary, NULL);
    */

    cl_int clStatus;

    ocl_scatter.lifetimes_kernel              = clCreateKernel(ocl_scatter.program_lifetimes, "lifetime", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create lifetime kernel.");

    ocl_integration.apply_max_lifetime = clCreateKernel(ocl_integration.program_integration, "apply_max_lifetime", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create apply_sigma_component kernel.");

    ocl_integration.apply_sigma_comp_kernel   = clCreateKernel(ocl_integration.program_integration, "apply_sigma_component", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create apply_sigma_component kernel.");

    ocl_integration.apply_simpson_weights     = clCreateKernel(ocl_integration.program_integration, "apply_simpson_weights", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create apply_simpson_weights kernel.");

    ocl_integration.integrate_position_kernel = clCreateKernel(ocl_integration.program_integration, "integrate_to_position", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create integrate_to_position kernel.");

    ocl_integration.sum_kernel                = clCreateKernel(ocl_integration.program_integration, "sum", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sum kernel.");

    ocl_scatter.particle_settings   = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ParticleSettings), nullptr, &clStatus);
    ocl_scatter.impurity_settings   = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ImpuritySettings), nullptr, &clStatus);
    ocl_scatter.simulation_settings = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(SimulationSettings), nullptr, &clStatus);
}

SimulationCL::~SimulationCL()
{
    clReleaseMemObject(ocl_scatter.particle_settings);
    clReleaseMemObject(ocl_scatter.impurity_settings);
    clReleaseMemObject(ocl_scatter.simulation_settings);
    clReleaseMemObject(ocl_scatter.impurities);
    clReleaseKernel(ocl_scatter.lifetimes_kernel);
    clReleaseProgram(ocl_scatter.program_lifetimes);

    clReleaseMemObject(ocl_integration.lifetimes_sigma_xx);
    clReleaseMemObject(ocl_integration.lifetimes_sigma_xy);
    clReleaseMemObject(ocl_integration.incomplete_sum);
    clReleaseKernel(ocl_integration.apply_sigma_comp_kernel);
    clReleaseKernel(ocl_integration.apply_simpson_weights);
    clReleaseKernel(ocl_integration.integrate_position_kernel);
    clReleaseKernel(ocl_integration.sum_kernel);
    clReleaseProgram(ocl_integration.program_integration);

    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}
