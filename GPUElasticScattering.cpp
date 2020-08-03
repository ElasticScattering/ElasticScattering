#include <windows.h>

#include <random>
#include <limits>
#include <vector>

#include "ElasticScattering.h"

// typedef unsigned int uint;

void GPUElasticScattering::Init(SimulationParameters sp) // todo arguments
{
    char		deviceStr[256];
    char		vendorStr[256];
    const char* source_file = "program.cl";

    double total_time;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    InitializeOpenCL(&ocl.deviceID, &ocl.context, &ocl.queue);

    PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    QueryPerformanceCounter(&beginClock);
    CompileOpenCLProgram(ocl.deviceID, ocl.context, source_file, &ocl.program);
    QueryPerformanceCounter(&endClock);

    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Time to build OpenCL Program: " << total_time * 1000 << " ms" << std::endl;

    PrepareOpenCLKernels(impurities, sp.particle_count);

    impurities.clear();
    impurities.resize(sp.impurity_count);

    std::cout << "Impurity region: " << -sp.particle_speed * sp.tau << ", " << sp.region_size + sp.particle_speed * sp.tau << std::endl;
    std::uniform_real_distribution<double> unif(-sp.particle_speed * sp.tau, sp.region_size + sp.particle_speed * sp.tau);
    std::random_device r;
    std::default_random_engine re(0);

    for (int i = 0; i < sp.impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    QueryPerformanceCounter(&beginClock);
    //Compute();
    QueryPerformanceCounter(&endClock);
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;
}


void GPUElasticScattering::PrepareOpenCLKernels(std::vector<cl_double2> impurities, int particle_count)
{
    cl_int clStatus;

    ocl.kernel = clCreateKernel(ocl.program, "scatter", &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create kernel.");

    ocl.impb = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(cl_double2) * impurities.size(), nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impb, CL_TRUE, 0, sizeof(cl_double2) * impurities.size(), impurities.data(), 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't enqueue buffer.");

    ocl.alive_buffer = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(bool) * particle_count, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create imp buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 0, sizeof(cl_mem), (void*)&ocl.impb);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");

    double imp_radius = 0.003;
    clStatus = clSetKernelArg(ocl.kernel, 1, sizeof(double), (void*)&imp_radius);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 2, sizeof(cl_mem), (void*)&ocl.alive_buffer);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::Compute()
{
    cl_int clStatus;
    size_t global_work_size = (size_t)sp.particle_count;
    size_t local_work_size = 20;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.kernel, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    //clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    /*result = new bool[size];
    memset(result, 0, sizeof(bool) * size);
    clEnqueueReadBuffer(ocl.queue, ocl.db, CL_TRUE, 0, sizeof(bool) * size, result, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Failed to read back result.");
    */
}

GPUElasticScattering::~GPUElasticScattering()
{
    clReleaseMemObject(ocl.impb);
    clReleaseMemObject(ocl.alive_buffer);
    clReleaseKernel(ocl.kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}