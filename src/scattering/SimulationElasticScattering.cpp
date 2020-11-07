/*#include "windows.h"
#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;

    cl_kernel scatter_kernel;
    cl_kernel add_integral_weights_kernel;
    cl_kernel sum_kernel;
    
    cl_mem parameters;
    cl_mem impurities;
    cl_mem imp_index;

    cl_mem xx_buffer;
    cl_mem xy_buffer;
    
    cl_mem xx_sum;
    cl_mem xy_sum;
} OCLSimResources;

OCLSimResources ocl;

bool SimulationElasticScattering::Compute(ScatteringParameters& p_sp, v2& result)
{
    PrepareCompute(p_sp);

    size_t global_work_size[2] = { (size_t)sp.dim, (size_t)sp.dim };
    size_t local_work_size[2] = { min(sp.dim, 256), 256 / min(sp.dim, 256) };

    cl_int clStatus;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.scatter_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel execution.");

    clStatus = clSetKernelArg(ocl.add_integral_weights_kernel, 0, sizeof(cl_mem), (void*)&ocl.xx_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.add_integral_weights_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

    clStatus = clSetKernelArg(ocl.add_integral_weights_kernel, 0, sizeof(cl_mem), (void*)&ocl.xy_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.add_integral_weights_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

    const size_t half_size = particle_count / 2;
    const size_t max_work_items = min(sp.dim, 256);

    clStatus = clSetKernelArg(ocl.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl.xx_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl.xx_sum);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 2, sizeof(double) * min(sp.dim, 256), nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_kernel, 1, nullptr, &half_size, &max_work_items, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start sum_lifetimes kernel execution.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl.xy_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl.xy_sum);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 2, sizeof(double) * min(sp.dim, 256), nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_kernel, 1, nullptr, &half_size, &max_work_items, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start sum_lifetimes kernel execution.");

    clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    std::vector<double> results;
    results.resize(half_size / max_work_items);
    clEnqueueReadBuffer(ocl.queue, ocl.xx_sum, CL_TRUE, 0, sizeof(double) * half_size / max_work_items, results.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

    result.x = ComputeResult(results);

    results.clear();
    results.resize(half_size / max_work_items);
    clEnqueueReadBuffer(ocl.queue, ocl.xy_sum, CL_TRUE, 0, sizeof(double) * half_size / max_work_items, results.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

    result.y = ComputeResult(results);
}

bool SimulationElasticScattering::PrepareCompute(ScatteringParameters& p_sp)
{
    cl_int clStatus;

    CompleteSimulationParameters(p_sp);

    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed = (sp.dim != p_sp.dim);

    sp = p_sp;

    if (first_run || impurities_changed) {
        grid.Generate(sp);

        ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(v2) * grid.impurities.size(), nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create impurities buffer.");

        clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * grid.impurities.size(), grid.impurities.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't write to impurities buffer.");

        ocl.imp_index = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(int) * grid.imp_index.size(), nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create impurity index buffer.");

        clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(int) * grid.imp_index.size(), grid.imp_index.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't write to impurities buffer.");
    }

    if (first_run) {
        ocl.parameters = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(ScatteringParameters), nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");
    }

    if (first_run || work_size_changed) {
        ocl.xx_buffer = clCreateBuffer(ocl.context, CL_MEM_WRITE_ONLY, sizeof(double) * particle_count, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

        ocl.xy_buffer = clCreateBuffer(ocl.context, CL_MEM_WRITE_ONLY, sizeof(double) * particle_count, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

        ocl.xx_sum = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * particle_count / 2, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create summation buffer.");

        ocl.xy_sum = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * particle_count / 2, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create summation buffer.");
    }

    if (first_run) {
        ocl.scatter_kernel = clCreateKernel(ocl.program, "scatter_sim", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

        ocl.add_integral_weights_kernel = clCreateKernel(ocl.program, "add_integral_weights_2d", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

        ocl.sum_kernel = clCreateKernel(ocl.program, "sum", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    }

    first_run = false;

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(ScatteringParameters), (void*)&sp, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 0, sizeof(cl_mem), (void*)&ocl.parameters);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 1, sizeof(cl_mem), (void*)&ocl.impurities);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 2, sizeof(cl_mem), (void*)&ocl.imp_index);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 3, sizeof(cl_mem), (void*)&ocl.xx_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 4, sizeof(cl_mem), (void*)&ocl.xy_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.add_integral_weights_kernel, 0, sizeof(cl_mem), (void*)&ocl.xx_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    return true;
}

bool SimulationElasticScattering::Compute(ScatteringParameters& p_sp, double& result)
{
    return true;
}

SimulationElasticScattering::SimulationElasticScattering(bool use_gpu, bool show_info) 
{
    InitializeOpenCL(use_gpu, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (show_info)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    std::cout << "\nBuilding OpenCL program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "scatter.cl", &ocl.program);
    std::cout << "Compilation succeeded! Starting simulation..." << std::endl;
}
SimulationElasticScattering::SimulationElasticScattering() : SimulationElasticScattering(true, false) { }

SimulationElasticScattering::SimulationElasticScattering(const InitParameters& init) : SimulationElasticScattering(init.use_gpu, !init.dont_show_info) {}

SimulationElasticScattering::~SimulationElasticScattering()
{
    clReleaseMemObject(ocl.parameters);
    clReleaseMemObject(ocl.impurities);
    clReleaseMemObject(ocl.xx_buffer);
    clReleaseMemObject(ocl.xy_buffer);
    clReleaseMemObject(ocl.xx_sum);
    clReleaseMemObject(ocl.xy_sum);
    clReleaseKernel(ocl.scatter_kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}
*/