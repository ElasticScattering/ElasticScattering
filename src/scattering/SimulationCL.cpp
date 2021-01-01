#include "windows.h"
#include "Simulation.h"
#include "src/utils/OpenCLUtils.h"

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;

    //cl_kernel quadrant_lifetimes_kernel;
    //cl_kernel quadrant_apply_sigma_comp_kernel;
    cl_kernel lifetimes_kernel;
    cl_kernel apply_sigma_comp_kernel;
    cl_kernel integrate_particle_kernel;
    cl_kernel simpson_weight_particle_kernel;
    cl_kernel sum_kernel;

    cl_mem particle_settings;
    cl_mem impurity_settings;
    cl_mem simulation_settings;

    cl_mem impurities;
    cl_mem cell_indices;
    cl_mem metrics;

    cl_mem lifetimes;
    cl_mem lifetimes_particle;

    cl_mem sigma_xx;
    cl_mem sigma_xx_particle;
    cl_mem sigma_xy;
    cl_mem sigma_xy_particle;

    cl_mem incomplete_sum;
} OCLSimResources;

OCLSimResources ocl;



Sigma SimulationCL::DeriveTemperature(const double temperature) const
{
    Sigma ir;

    return ir;
}

IterationResult SimulationCL::DeriveTemperatureWithImages(const double temperature) const
{

    IterationResult ir;

    return ir;

    /*
    cl_int clStatus;

    // temperature is niet genoeg, moet sp meegeven...
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(ScatteringParameters), (void*)&sp, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    // Set up kernels
    clStatus = clSetKernelArg(ocl.quadrant_apply_sigma_comp_kernel, 0, sizeof(cl_mem), (void*)&ocl.lifetimes);
    clStatus = clSetKernelArg(ocl.quadrant_apply_sigma_comp_kernel, 1, sizeof(cl_mem), (void*)&ocl.parameters);

    // Integrate particle kernel.
    // First argument should be set before execution.
    clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 1, sizeof(cl_mem), (void*)&ocl.parameters);
    clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 2, sizeof(cl_mem), (void*)&ocl.lifetimes_particle);

    clStatus = clSetKernelArg(ocl.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl.incomplete_sum);
    clStatus = clSetKernelArg(ocl.sum_kernel, 2, sizeof(double) * items_in_workgroup, nullptr);

    int mode = 0;
    clStatus = clSetKernelArg(ocl.quadrant_apply_sigma_comp_kernel, 2, sizeof(int), (void*)&mode);
    clStatus = clSetKernelArg(ocl.quadrant_apply_sigma_comp_kernel, 3, sizeof(cl_mem), (void*)&ocl.sigma_xx);

    mode = 1;
    clStatus = clSetKernelArg(ocl.quadrant_apply_sigma_comp_kernel, 2, sizeof(int), (void*)&mode);
    clStatus = clSetKernelArg(ocl.quadrant_apply_sigma_comp_kernel, 3, sizeof(cl_mem), (void*)&ocl.sigma_xy);


    ocl.sigma_xx                       = clCreateBuffer(ocl.context, CL_MEM_WRITE_ONLY, sizeof(double) * particle_count, nullptr, &clStatus);
    ocl.sigma_xy                       = clCreateBuffer(ocl.context, CL_MEM_WRITE_ONLY, sizeof(double) * particle_count, nullptr, &clStatus);
    ocl.incomplete_sum                 = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * particle_count / 2, nullptr, &clStatus);

    const size_t items_in_work_group = min(sp.dim, 256);

    PrepareKernels(sp, items_in_work_group);

    cl_int clStatus;

    const size_t particles_work_size = sp.dim * sp.dim;
    const size_t total_work_size = particles_work_size * sp.values_per_particle;
    const size_t half_size = total_work_size / 2;
    size_t global_work_size[3] = { (size_t)sp.dim, (size_t)sp.dim, 4 };
    size_t local_work_size[3] = { items_in_work_group, 256 / items_in_work_group, 1 };

    std::vector<double> particle_lifetimes(particles_work_size);
    {
        //	const double integrand_factor = sp->integrand_angle_area / ((values_per_quadrant - 1) * 3.0);

        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 0, sizeof(cl_mem), (void*)&ocl.lifetimes);
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.integrate_particle_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start integration kernel execution.");

        clEnqueueReadBuffer(ocl.queue, ocl.lifetimes_particle, CL_TRUE, 0, sizeof(double) * total_work_size, particle_lifetimes.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
    }

    std::vector<double> sigma_xx(particles_work_size);
    std::vector<double> sigma_xy(particles_work_size);
    {
        // Sigma lifetimes
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.quadrant_apply_sigma_comp_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel execution.");

        // Sigma xx per particle
        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 0, sizeof(cl_mem), (void*)&ocl.sigma_xx);
        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 2, sizeof(cl_mem), (void*)&ocl.sigma_xx_particle);
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.integrate_particle_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start integration kernel execution.");

        // Sigma xy per particle
        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 0, sizeof(cl_mem), (void*)&ocl.sigma_xy);
        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 2, sizeof(cl_mem), (void*)&ocl.sigma_xy_particle);
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.integrate_particle_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start integration kernel execution.");

        clEnqueueReadBuffer(ocl.queue, ocl.sigma_xx_particle, CL_TRUE, 0, sizeof(double) * total_work_size, sigma_xx.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back sigma xx.");

        clEnqueueReadBuffer(ocl.queue, ocl.sigma_xy_particle, CL_TRUE, 0, sizeof(double) * total_work_size, sigma_xy.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back sigma xy.");
    }

    Sigma sigma;
    {
        // Integrate sigmas
        // Sigma XX
        std::vector<double> incomplete_sum;
        incomplete_sum.resize(half_size / items_in_work_group);

        clStatus = clSetKernelArg(ocl.simpson_weight_particle_kernel, 0, sizeof(cl_mem), (void*)&ocl.sigma_xx_particle);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.simpson_weight_particle_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

        clStatus = clSetKernelArg(ocl.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl.sigma_xx_particle);
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);

        clEnqueueReadBuffer(ocl.queue, ocl.incomplete_sum, CL_TRUE, 0, sizeof(double) * half_size / items_in_work_group, incomplete_sum.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

        for (int i = 0; i < incomplete_sum.size(); i++)
            sigma.xx += incomplete_sum[i];
    }

    {
        // Sigma XY
        std::vector<double> incomplete_sum;
        incomplete_sum.resize(half_size / items_in_work_group);

        clStatus = clSetKernelArg(ocl.simpson_weight_particle_kernel, 0, sizeof(cl_mem), (void*)&ocl.sigma_xy_particle);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.simpson_weight_particle_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

        clStatus = clSetKernelArg(ocl.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl.sigma_xy_particle);
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);

        clEnqueueReadBuffer(ocl.queue, ocl.incomplete_sum, CL_TRUE, 0, sizeof(double) * half_size / items_in_work_group, incomplete_sum.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

        for (int i = 0; i < incomplete_sum.size(); i++)
            sigma.xy += incomplete_sum[i];
    }
    */
}

void SimulationCL::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
    cl_int clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.particle_settings, CL_TRUE, 0, sizeof(ParticleSettings), (void*)&ps, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurity_settings, CL_TRUE, 0, sizeof(ImpuritySettings), (void*)&is, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.simulation_settings, CL_TRUE, 0, sizeof(SimulationSettings), (void*)&ss, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.metrics, CL_TRUE, 0, sizeof(Metrics), (void*)&metrics, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    {
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 0, sizeof(cl_mem), (void*)&ocl.simulation_settings);
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 1, sizeof(cl_mem), (void*)&ocl.particle_settings);
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 2, sizeof(cl_mem), (void*)&ocl.impurity_settings);
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 3, sizeof(cl_mem), (void*)&ocl.impurities);
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 4, sizeof(cl_mem), (void*)&ocl.cell_indices);
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 5, sizeof(cl_mem), (void*)&ocl.lifetimes);
        clStatus = clSetKernelArg(ocl.lifetimes_kernel, 6, sizeof(cl_mem), (void*)&ocl.metrics);
    }
    
    size_t global_work_size[3] = { (size_t)ss.particles_per_row, (size_t)ss.particles_per_row, ss.values_per_particle };
    const size_t values_in_work_group = min(ss.values_per_particle, 256);
    size_t local_work_size[3] = { 256 / values_in_work_group, 1, values_in_work_group }; // @Todo, dit verifiëren!

    //const size_t items_in_work_group = min(ss.particles_per_row, 256);
    //size_t local_work_size[3] = { items_in_work_group, 256 / items_in_work_group, 1 };

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.lifetimes_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel.");
    clFinish(ocl.queue);

    std::vector<Metrics> metrics_holder(1);
    clEnqueueReadBuffer(ocl.queue, ocl.metrics, CL_TRUE, 0, sizeof(Metrics), metrics_holder.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
    metrics = metrics_holder[0];
}

void SimulationCL::UploadImpurities(const Grid& grid)
{
    cl_int clStatus;

    ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(v2) * grid.GetImpurities().size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create impurities buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * grid.GetImpurities().size(), grid.GetImpurities().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't upload impurities.");

    ocl.cell_indices = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(int) * grid.GetIndex().size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create cell index buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.cell_indices, CL_TRUE, 0, sizeof(int) * grid.GetIndex().size(), grid.GetIndex().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't upload cell indices.");
}

void SimulationCL::PrepareKernels(const Settings& ss, const size_t items_in_workgroup)
{
}

SimulationCL::SimulationCL() : SimulationCL(true, true) {};

SimulationCL::SimulationCL(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant)
{
    InitializeOpenCL(true, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (true)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    std::cout << "\nBuilding OpenCL program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "src/scattering/cl/scatter.cl", &ocl.program);
    std::cout << "Compilation succeeded! Starting simulation..." << std::endl;

    cl_int clStatus;

    ocl.lifetimes_kernel      = clCreateKernel(ocl.program, "lifetime", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    /*
    ocl.quadrant_apply_sigma_comp_kernel         = clCreateKernel(ocl.program, "quadrant_sigma_lifetimes", &clStatus);
    ocl.integrate_particle_kernel      = clCreateKernel(ocl.program, "integrate_particle", &clStatus);
    ocl.simpson_weight_particle_kernel = clCreateKernel(ocl.program, "simpson_weight_particle_kernel", &clStatus);
    ocl.sum_kernel                     = clCreateKernel(ocl.program, "sum", &clStatus);
    */
 
    ocl.particle_settings   = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ParticleSettings), nullptr, &clStatus);
    ocl.impurity_settings   = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ImpuritySettings), nullptr, &clStatus);
    ocl.simulation_settings = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(SimulationSettings), nullptr, &clStatus);
}

SimulationCL::~SimulationCL()
{
    clReleaseMemObject(ocl.particle_settings);
    clReleaseMemObject(ocl.impurity_settings);
    clReleaseMemObject(ocl.simulation_settings);
    clReleaseMemObject(ocl.impurities);
    clReleaseMemObject(ocl.sigma_xx);
    clReleaseMemObject(ocl.sigma_xy);
    clReleaseMemObject(ocl.incomplete_sum);
    
    clReleaseKernel(ocl.lifetimes_kernel);
    clReleaseKernel(ocl.apply_sigma_comp_kernel);
    clReleaseKernel(ocl.integrate_particle_kernel);
    clReleaseKernel(ocl.simpson_weight_particle_kernel);
    clReleaseKernel(ocl.sum_kernel);
    
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}
