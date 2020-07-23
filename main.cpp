#include <math.h>
#include <windows.h>
#include <fstream>
#include <sstream>

#include <random>
#include <limits>
#include <vector>
#include <unordered_map>

#include "opencl_utils.h"

typedef unsigned int uint;

#define RUN_CPU_SIM

enum class Mode {
    LIFETIME,
    STATS,
    AVG_DISTANCE,
    CONDUCTIVITY
}; 

typedef struct
{
    int num_iterations;
    Mode mode;
    bool show_info;
} InitParameters;

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

cl_int clStatus;
double *result; // @Temporary

void GPUElasticScattering(OCLResources *p_ocl, size_t size)
{
    size_t local_work_size = 20;

    clStatus = clEnqueueNDRangeKernel(p_ocl->queue, p_ocl->kernel, 1, nullptr, &size, &local_work_size, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    //clStatus = clFinish(p_ocl->queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    result = new double[size];
    memset(result, 0, sizeof(double) * size);
    clEnqueueReadBuffer(p_ocl->queue, p_ocl->db, CL_TRUE, 0, sizeof(double) * size, result, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Failed to read back result.");
}

void ParseArgs(OCLResources* p_ocl, int argc, char** argv, InitParameters *p_init) {
    if (argc != 4) {
        std::cout << "Usage: ElasticScattering [num iterations] [lifetime | distance | stats | conductivity] [show | no-show]" << std::endl;
        exit(0);
    }

    p_init->num_iterations = atoi(argv[1]);

    std::unordered_map<std::string, Mode> modes
    {
        {"lifetime",     Mode::LIFETIME},
        {"distance",     Mode::AVG_DISTANCE},
        {"stats",        Mode::STATS},
        {"conductivity", Mode::CONDUCTIVITY}
    };
    std::string key;
    key.assign(argv[2], strlen(argv[2]));
    auto iterator = modes.find(key);
    ERR_FAIL_COND_MSG(iterator == modes.end(), "Couldn't understand second command line argument.");
    p_init->mode = iterator->second;
    
    p_init->show_info = strcmp(argv[3], "show");
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

    std::cout << "+---------------------------------------------------+" << std::endl;
}

int main(int argc, char* argv[]) 
{
    OCLResources ocl;
    char		deviceStr[256];
    char		vendorStr[256];
    const char* source_file = "program.cl";

    uint startHeight = 32, startWidth = 32;
    uint particleCount = 1'000'000;
    uint impurityCount = 1000;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    InitParameters init;
    ParseArgs(&ocl, argc, argv, &init);

    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
    std::cout << "Initial field size: " << startHeight << ", " << startWidth << std::endl;

    // Initialize buffers.
    std::uniform_real_distribution<double> unif(0, 1000);
    std::default_random_engine re;
    int size = particleCount;
    double* data = new double[size];
    for (int i = 0; i < size; i++)
        data[i] = pow(10, 200); // unif(re);

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
    
    QueryPerformanceCounter(&beginClock);
    GPUElasticScattering(&ocl, particleCount);
    QueryPerformanceCounter(&endClock);
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;

    std::cout.precision(64); // std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < min(size -1, 100); i++)
        std::cout << result[i] << ", ";

    std::cout << result[size - 1] << std::endl;

    Cleanup(&ocl);

    return 0;
}
