#include <windows.h>

#include <random>
#include <limits>
#include <vector>
#include <unordered_map>

#include "ElasticScattering.h"

typedef unsigned int uint;

#define RUN_CPU_SIM

cl_int clStatus;
bool* result; // @Temporary

void ElasticScattering::GPUElasticScattering(size_t size)
{
    size_t local_work_size = 20;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.kernel, 1, nullptr, &size, &local_work_size, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    //clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    result = new bool[size];
    memset(result, 0, sizeof(bool) * size);
    clEnqueueReadBuffer(ocl.queue, ocl.db, CL_TRUE, 0, sizeof(bool) * size, result, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Failed to read back result.");
}

void ElasticScattering::ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->num_iterations = 1;
    p_init->mode = Mode::LIFETIME;
    p_init->show_info = true;
    return;

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

void ElasticScattering::PrepareOpenCLKernels()
{
    ocl.kernel = clCreateKernel(ocl.program, "scatter", &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create kernel.");

    ocl.impb = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(cl_double2) * impurity_count, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impb, CL_TRUE, 0, sizeof(cl_double2) * impurity_count, imp_data, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't enqueue buffer.");

    ocl.alive_buffer = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(bool) * particle_count, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.alive_buffer, CL_TRUE, 0, sizeof(bool) * particle_count, alive_data, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't enqueue buffer.");


    clStatus = clSetKernelArg(ocl.kernel, 0, sizeof(cl_mem), (void*)&ocl.impb);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");

    double imp_radius = 0.003;
    clStatus = clSetKernelArg(ocl.kernel, 1, sizeof(double), (void*)&imp_radius);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.kernel, 2, sizeof(cl_mem), (void*)&ocl.alive_buffer);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't set argument to buffer.");

}

void ElasticScattering::Cleanup()
{
    clReleaseMemObject(ocl.impb);
    clReleaseMemObject(ocl.alive_buffer);
    clReleaseKernel(ocl.kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}

void ElasticScattering::Init(int argc, char* argv[])
{
    char		deviceStr[256];
    char		vendorStr[256];
    const char* source_file = "program.cl";

    uint startHeight = 32, startWidth = 32;
    particle_count = 1'000'000;
    impurity_count = 1000;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    InitParameters init;
    ParseArgs(argc, argv, &init);

    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
    std::cout << "Initial field size: " << startHeight << ", " << startWidth << std::endl;

    // Initialize buffers.
    std::uniform_real_distribution<double> unif(0, 1000);
    std::default_random_engine re;
    imp_data = new cl_double2[impurity_count];
    for (int i = 0; i < impurity_count; i++)
    {
        imp_data[i].x = unif(re);
        imp_data[i].y = unif(re);
    }

    alive_data = (bool*)malloc(particle_count * sizeof(bool));
    ERR_FAIL_COND_MSG(!alive_data, "Could not init arrays.")
    memset(alive_data, false, particle_count * sizeof(bool));

#ifdef RUN_CPU_SIM
#endif

    // Setup
    InitializeOpenCL(&ocl.deviceID, &ocl.context, &ocl.queue);

    PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    QueryPerformanceCounter(&beginClock);
    CompileOpenCLProgram(ocl.deviceID, ocl.context, source_file, &ocl.program);
    QueryPerformanceCounter(&endClock);

    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Time to build OpenCL Program: " << total_time * 1000 << " ms" << std::endl;

    PrepareOpenCLKernels();

    QueryPerformanceCounter(&beginClock);
    GPUElasticScattering(particle_count);
    QueryPerformanceCounter(&endClock);
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;

    std::cout.precision(64); // std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < min(particle_count - 1, 100); i++)
        std::cout << "(" << result[i] << "), ";

    //std::cout << "(" << result[particle_count - 1].x << ", " << result[particle_count - 1].y << + ")" << std::endl;
}