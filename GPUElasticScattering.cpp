#include <windows.h>

#include <random>
#include <limits>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"
#include "Details.h"

// typedef unsigned int uint;
typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;
    cl_kernel kernel;

    cl_mem db;
    cl_mem impb;
    cl_mem lifetimes;
} OCLResources;

OCLResources ocl;


void GPUElasticScattering::Init(SimulationParameters p_sp) // todo arguments
{
    char		deviceStr[256];
    char		vendorStr[256];
    const char* source_file = "scatterB.cl";
    sp = p_sp;
    if (sp.angular_speed == 0) sp.particle_max_lifetime = sp.tau;
    else {
        double bound_time = GetBoundTime(sp.phi, sp.alpha, sp.angular_speed, true, false);
        sp.particle_max_lifetime = MIN(sp.tau, bound_time);
    }
    std::cout << "Particle max lifetime: " << sp.particle_max_lifetime << std::endl;

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


    impurities.clear();
    impurities.resize(sp.impurity_count);

    std::cout << "Impurity region: " << -sp.particle_speed * sp.tau << ", " << sp.region_size + sp.particle_speed * sp.tau << std::endl;
    std::uniform_real_distribution<double> unif(-sp.particle_speed * sp.tau, sp.region_size + sp.particle_speed * sp.tau);
    std::random_device r;
    std::default_random_engine re(0);

    for (int i = 0; i < sp.impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    sp.particle_count *= 100;
    sp.particle_row_count = sqrt(sp.particle_count);
//    lifetimes.clear();
  //  lifetimes.resize(sp.particle_count, 0);

    PrepareOpenCLKernels();
}


void GPUElasticScattering::PrepareOpenCLKernels()
{
    cl_int clStatus;

    ocl.kernel = clCreateKernel(ocl.program, "scatterB", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.impb = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(v2) * impurities.size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impb, CL_TRUE, 0, sizeof(v2) * impurities.size(), impurities.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

    ocl.lifetimes = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * sp.particle_count, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

    
    
    clStatus = clSetKernelArg(ocl.kernel, 0, sizeof(double), (void*)&sp.region_size);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 1, sizeof(double), (void*)&sp.particle_max_lifetime);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 2, sizeof(double), (void*)&sp.particle_speed);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 3, sizeof(double), (void*)&sp.particle_mass);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 4, sizeof(double), (void*)&sp.impurity_radius);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 5, sizeof(double), (void*)&sp.tau);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 6, sizeof(double), (void*)&sp.alpha);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 7, sizeof(double), (void*)&sp.phi);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 8, sizeof(double), (void*)&sp.magnetic_field);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 9, sizeof(double), (void*)&sp.angular_speed);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 10, sizeof(int), (void*)&sp.impurity_count);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 11, sizeof(cl_mem), (void*)&ocl.impb);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
    
    clStatus = clSetKernelArg(ocl.kernel, 12, sizeof(cl_mem), (void*)&ocl.lifetimes);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::Compute()
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    QueryPerformanceCounter(&beginClock);

    cl_int clStatus;
    size_t global_work_size = (size_t)sp.particle_count;
    size_t local_work_size = 20;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.kernel, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    
    lifetimes.clear();
    lifetimes.resize(sp.particle_count, 0);

    clEnqueueReadBuffer(ocl.queue, ocl.lifetimes, CL_TRUE, 0, sizeof(double) * sp.particle_count, lifetimes.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus != CL_SUCCESS, clStatus, "Failed to read back result.");

    QueryPerformanceCounter(&endClock);
    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;


    std::cout << "\n\Results:" << std::endl;
    for (int i = 0; i < MIN(lifetimes.size(), 200); i++)
        std::cout << lifetimes[i] << ", ";
    std::cout << "..." << std::endl;

    MakeTexture(sp);
}

GPUElasticScattering::~GPUElasticScattering()
{
    clReleaseMemObject(ocl.impb);
    clReleaseMemObject(ocl.lifetimes);
    clReleaseKernel(ocl.kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}