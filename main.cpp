#include <math.h>
#include <windows.h>
#include <fstream>
#include <sstream>

#include <random>
#include <limits>
#include <vector>

#include "opencl_utils.h"

typedef unsigned int uint;

#define RUN_CPU_SIM

typedef struct
{
    // CL platform handles:
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;
    cl_kernel kernel;
    
    // CL memory buffers
    cl_mem db;
    cl_mem impb;
} OCLResources;

cl_platform_id selected_platform;
cl_device_id selected_device;

cl_int clStatus;
double *result;

void GPUElasticScattering(OCLResources *p_ocl, size_t size)
{
    size_t local_work_size = 20;

    clStatus = clEnqueueNDRangeKernel(p_ocl->queue, p_ocl->kernel, 1, nullptr, &size, &local_work_size, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    clStatus = clFinish(p_ocl->queue);

    result = new double[size];
    memset(result, 0, sizeof(double) * size);
    clEnqueueReadBuffer(p_ocl->queue, p_ocl->db, CL_TRUE, 0, sizeof(double) * size, result, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Failed to read back result.");
}

void ParseArgs(OCLResources* p_ocl, int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: ElasticScattering [lifetime | distance | stats | conductivity] [show | no-show]" << std::endl;
        exit(0);
    }
}

void PrepareOpenCLKernels(OCLResources* p_ocl, size_t size, double *data)
{
    p_ocl->kernel = clCreateKernel(p_ocl->program, "double_precision", &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create kernel.");

    p_ocl->db = clCreateBuffer(p_ocl->context, CL_MEM_READ_WRITE, sizeof(double) * size, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create buffer.");

    clStatus = clEnqueueWriteBuffer(p_ocl->queue, p_ocl->db, CL_TRUE, 0, sizeof(double) * size, data, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't enqueue buffer.");
    
    clStatus = clSetKernelArg(p_ocl->kernel, 0, sizeof(cl_mem), (void *)&p_ocl->db);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");
}

void Cleanup(OCLResources* p_ocl)
{
    clReleaseMemObject(p_ocl->db);
    clReleaseKernel(p_ocl->kernel);
    clReleaseProgram(p_ocl->program);
    clReleaseCommandQueue(p_ocl->queue);
    clReleaseContext(p_ocl->context);

    std::cout << "-----------------------------------------------" << std::endl;
}

int main(int argc, char* argv[]) 
{
    OCLResources ocl;
    char		deviceStr[256];
    char		vendorStr[256];
    const char* source_file = "program.cl";

    uint startHeight = 32, startWidth = 32;
    uint impurityCount = 100;
    uint particleCount = 10000; //100'000'000; // batch size
    uint impurityCount = 1000;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    //parse_args(ocl, argc, argv);

    std::cout << "\n\n-----------------------------------------------" << std::endl;
    std::cout << "Initial field size: " << startHeight << ", " << startWidth << std::endl;

    // Initialize buffers.
    std::uniform_real_distribution<double> unif(0, 1000);
    std::default_random_engine re;
    int size = particleCount;
    double* data = new double[size];
    for (int i = 0; i < size; i++)
        data[i] = unif(re);

#ifdef RUN_CPU_SIM
#endif
    
    // Setup
    InitializeOpenCL(deviceStr, vendorStr, &ocl.deviceID, &ocl.context, &ocl.queue);

    PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    QueryPerformanceCounter(&beginClock);
    CompileOpenCLProgram(ocl.deviceID, ocl.context, source_file, &ocl.program);
    QueryPerformanceCounter(&endClock);

    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Time to build OpenCL Program: " << total_time * 1000 << " ms" << std::endl;

    PrepareOpenCLKernels(&ocl, particleCount, data);
    
    // Launch kernel
    QueryPerformanceCounter(&beginClock);
    GPUElasticScattering(&ocl, particleCount);
    QueryPerformanceCounter(&endClock);
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;

    /*
    std::cout.precision(64); // std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < size - 1; i++)
        std::cout << result[i] << ", ";

    std::cout << result[size - 1] << std::endl;
    */

    Cleanup(&ocl);

    return 0;
}
