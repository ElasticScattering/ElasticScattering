#include "windows.h"
#include "ElasticScattering.h"
#include "src/utils/OpenCLUtils.h"

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;

    cl_kernel quadrant_lifetimes_kernel;
    cl_kernel quadrant_sigmas_kernel;
    cl_kernel integrate_particle_kernel;
    cl_kernel simpson_weight_particle_kernel;
    cl_kernel sum_kernel;

    cl_mem parameters;
    cl_mem impurities;
    cl_mem imp_index;

    cl_mem lifetimes;
    cl_mem lifetimes_particle;

    cl_mem sigma_xx;
    cl_mem sigma_xx_particle;
    cl_mem sigma_xy;
    cl_mem sigma_xy_particle;

    cl_mem incomplete_sum;
} OCLSimResources;

OCLSimResources ocl;

IterationResult ElasticScatteringCL::ComputeIteration(const ScatteringParameters& sp, const ImpurityIndex& grid)
{
    const size_t items_in_work_group = min(sp.dim, 256);

    UploadImpurities(grid);
    PrepareKernels(sp, items_in_work_group);

    cl_int clStatus;

    const size_t particles_work_size = sp.dim * sp.dim;
    const size_t total_work_size = particles_work_size * sp.values_per_particle;
    const size_t half_size = total_work_size / 2;
    size_t global_work_size[3] = { (size_t)sp.dim, (size_t)sp.dim, 4 };
    size_t local_work_size[3] = { items_in_work_group, 256 / items_in_work_group, 1 };

    std::vector<double> particle_lifetimes(particles_work_size);
    {
        // Lifetimes
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.quadrant_lifetimes_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel execution.");

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
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.quadrant_sigmas_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
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

    IterationResult ir;
    ir.particle_lifetimes = particle_lifetimes;
    ir.sigmas.xx_buffer = sigma_xx;
    ir.sigmas.xx_buffer = sigma_xy;
    ir.result = sigma;

    return ir;
}

void ElasticScatteringCL::UploadImpurities(const ImpurityIndex& grid)
{
    cl_int clStatus;

    ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(v2) * grid.impurity_count, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create impurities buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * grid.impurity_count, grid.GetImpurities().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't write to impurities buffer.");

    ocl.imp_index = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(int) * grid.GetIndex().size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create impurity index buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(int) * grid.GetIndex().size(), grid.GetIndex().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't write to impurities buffer.");
}

void ElasticScatteringCL::PrepareKernels(const ScatteringParameters& sp, const size_t items_in_workgroup)
{
    cl_int clStatus;

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(ScatteringParameters), (void*)&sp, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    {
        // Scatter kernel.
        clStatus = clSetKernelArg(ocl.quadrant_lifetimes_kernel, 0, sizeof(cl_mem), (void*)&ocl.parameters);
        clStatus = clSetKernelArg(ocl.quadrant_lifetimes_kernel, 1, sizeof(cl_mem), (void*)&ocl.impurities);
        clStatus = clSetKernelArg(ocl.quadrant_lifetimes_kernel, 2, sizeof(cl_mem), (void*)&ocl.imp_index);
        clStatus = clSetKernelArg(ocl.quadrant_lifetimes_kernel, 3, sizeof(cl_mem), (void*)&ocl.lifetimes);
    }

    {
        // Sigma kernel.
        clStatus = clSetKernelArg(ocl.quadrant_sigmas_kernel, 0, sizeof(cl_mem), (void*)&ocl.lifetimes);
        clStatus = clSetKernelArg(ocl.quadrant_sigmas_kernel, 1, sizeof(cl_mem), (void*)&ocl.parameters);
        clStatus = clSetKernelArg(ocl.quadrant_sigmas_kernel, 2, sizeof(cl_mem), (void*)&ocl.sigma_xx);
        clStatus = clSetKernelArg(ocl.quadrant_sigmas_kernel, 3, sizeof(cl_mem), (void*)&ocl.sigma_xy);
    }

    {
        // Integrate particle kernel.
        // First argument should be set before execution.
        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 1, sizeof(cl_mem), (void*)&ocl.parameters);
        clStatus = clSetKernelArg(ocl.integrate_particle_kernel, 2, sizeof(cl_mem), (void*)&ocl.lifetimes_particle);
    }

    {
        clStatus = clSetKernelArg(ocl.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl.incomplete_sum);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.sum_kernel, 2, sizeof(double) * items_in_workgroup, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
    }
}

ElasticScatteringCL::ElasticScatteringCL(bool use_gpu, bool show_info, int particle_count)
{
    InitializeOpenCL(use_gpu, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (show_info)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    std::cout << "\nBuilding OpenCL program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "scatter.cl", &ocl.program);
    std::cout << "Compilation succeeded! Starting simulation..." << std::endl;


    cl_int clStatus;

    ocl.quadrant_lifetimes_kernel = clCreateKernel(ocl.program, "quadrant_lifetime", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.quadrant_sigmas_kernel = clCreateKernel(ocl.program, "quadrant_sigma_lifetimes", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.integrate_particle_kernel = clCreateKernel(ocl.program, "integrate_particle", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.simpson_weight_particle_kernel = clCreateKernel(ocl.program, "simpson_weight_particle_kernel", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.sum_kernel = clCreateKernel(ocl.program, "sum", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.parameters = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ScatteringParameters), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

    ocl.sigma_xx = clCreateBuffer(ocl.context, CL_MEM_WRITE_ONLY, sizeof(double) * particle_count, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create particle_lifetimes buffer.");

    ocl.sigma_xy = clCreateBuffer(ocl.context, CL_MEM_WRITE_ONLY, sizeof(double) * particle_count, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create particle_lifetimes buffer.");

    ocl.incomplete_sum = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * particle_count / 2, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sum buffer.");
}

ElasticScatteringCL::~ElasticScatteringCL()
{
    clReleaseMemObject(ocl.parameters);
    clReleaseMemObject(ocl.impurities);
    clReleaseMemObject(ocl.sigma_xx);
    clReleaseMemObject(ocl.sigma_xy);
    clReleaseMemObject(ocl.incomplete_sum);
    clReleaseKernel(ocl.quadrant_lifetimes_kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}
