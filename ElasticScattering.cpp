#include <windows.h>

#include <random>
#include <limits>
#include <vector>
#include <unordered_map>

#include "ElasticScattering.h"

typedef unsigned int uint;

const double PI = 3.141592653589793238463;
#define RUN_CPU_SIM


void ElasticScattering::CPUElasticScattering2(const SimulationParameters sp, const cl_double2* imp_pos, cl_double* lifetime_results)
{
    const int particles_in_row = sqrt(sp.particle_count);

    const bool clockwise = true;

    // Map phi to the interval[-alpha, 2pi - alpha).
    double phi2 = fmod(sp.phi + sp.alpha, (PI * 2.0)) - sp.alpha;

    // Map to the lower bound, so -alpha + n pi/2
    double low_bound = floor((phi2 + sp.alpha) / (PI * 0.5)) * (PI * 0.5);

    // Remaining is the distance to this boundary in rad.
    double remaining = phi2 - low_bound + sp.alpha;

    double dphi;
    if (!clockwise && false) dphi = remaining;
    else if (!clockwise)     dphi = 2*sp.alpha - remaining;
    else if (!false)         dphi = 2 * sp.alpha - remaining;
    else                     dphi = remaining;
    
    double bound_time = max(sp.tau, dphi / sp.angular_speed);
    
    for (int j = 0; j < particles_in_row; j++)
    {
        for (int i = 0; i < particles_in_row; i++)
        {
            cl_double2 pos;
            pos.x = sp.region_size * (double(i) / particles_in_row);
            pos.y = sp.region_size * (double(j) / particles_in_row);

            cl_double2 vel = { sp.particle_speed * cos(sp.phi), sp.particle_speed * sin(sp.phi) };

            double lifetime = bound_time;
        }
    }
}

void ElasticScattering::CPUElasticScattering(const SimulationParameters sp, const cl_double2 *imp_pos, cl_double *lifetime_results) 
{
    const int particles_in_row = sqrt(sp.particle_count);

    for (int j = 0; j < particles_in_row; j++) 
    {
        for (int i = 0; i < particles_in_row; i++)
        {
            cl_double2 pos;
            pos.x = sp.region_size * (double(i) / particles_in_row);
            pos.y = sp.region_size * (double(j) / particles_in_row);

            cl_double2 vel = { sp.particle_speed * cos(sp.phi), sp.particle_speed * sin(sp.phi) };
            
            double lifetime = sp.tau;

            for (int k = 0; k < impurity_count; k++)
            {
                const cl_double2 ip = imp_pos[k];
                const cl_double2 unit = { vel.x / sp.particle_speed, vel.y / sp.particle_speed };
                const cl_double2 projected = { pos.x + (ip.x - pos.x) * unit.x, pos.y + (ip.y - pos.y) * unit.y };

                const double a = pow(projected.x - ip.x, 2.0) + pow(projected.y - ip.y, 2.0);
                if (a > sp.impurity_radius_sq) {
                    continue; //@Speedup, if distance is greater than current min continue as well.
                }

                double L = sqrt(sp.impurity_radius_sq - a);

                cl_double2 time_taken;
                if (vel.x != 0) time_taken = { ((projected.x - L * unit.x) - pos.x) / vel.x, ((projected.x + L * unit.x) - pos.x) / vel.x };
                else            time_taken = { ((projected.y - L * unit.y) - pos.y) / vel.y, ((projected.y + L * unit.y) - pos.y) / vel.y };


                if ((time_taken.s0 * time_taken.s1) < 0) {
                    lifetime = 0;
                    break;
                }

                if (time_taken.s0 > 0 && time_taken.s0 < lifetime) {
                    lifetime = time_taken.s0;
                }
                if (time_taken.s1 > 0 && time_taken.s1 < lifetime) {
                    lifetime = time_taken.s1;
                }
            }

            lifetime_results[j * particles_in_row + i] = lifetime;
        }
    }
}

void ElasticScattering::GPUElasticScattering(size_t size)
{
    size_t local_work_size = 20;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.kernel, 1, nullptr, &size, &local_work_size, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    //clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    /*result = new bool[size];
    memset(result, 0, sizeof(bool) * size);
    clEnqueueReadBuffer(ocl.queue, ocl.db, CL_TRUE, 0, sizeof(bool) * size, result, 0, nullptr, nullptr);
    CL_ERR_FAIL_COND_MSG(clStatus != CL_SUCCESS, clStatus, "Failed to read back result.");
    */
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

    double total_time;

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    InitParameters init;
    ParseArgs(argc, argv, &init);

    SimulationParameters sp;
    sp.region_size        = 5e-6;
    sp.particle_count     = 10'000; //100'000'000;
    sp.particle_speed     = 7e5;
    sp.impurity_count     = 1000;
    sp.impurity_radius    = 1.5e-8;
    sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
    sp.tau                = 1e-12;
    sp.particle_mass      = 5 * 9.1e-31; 
    sp.alpha              = 3.14159 / 4.0;
    sp.phi                = sp.alpha;
    sp.magnetic_field     = 0.3;
    sp.angular_speed      = 1.602e-19 * sp.magnetic_field / sp.particle_mass;

    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "Start region size: (" << sp.region_size << ", " << sp.region_size << ")" << std::endl;
    std::cout << "Particles:         " << sp.particle_count << std::endl;
    std::cout << "Particle speed:    " << sp.particle_speed << std::endl;
    std::cout << "Particle mass:     " << sp.particle_mass << std::endl;
    std::cout << "Impurities:        " << sp.impurity_count << std::endl;
    std::cout << "Impurity radius:   " << sp.impurity_radius << std::endl;
    std::cout << "Tau:               " << sp.tau << std::endl;
    std::cout << "Alpha:             " << sp.alpha << std::endl;
    std::cout << "Phi:               " << sp.phi << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    // Initialize buffers.
    std::uniform_real_distribution<double> unif(0, 5e-6);
    std::random_device r;
    std::default_random_engine re(r());
    imp_data = new cl_double2[impurity_count];
    for (int i = 0; i < impurity_count; i++)
        imp_data[i] = { unif(re), unif(re) };

    double *lifetime_results = (double*)malloc(sp.particle_count * sizeof(double));
    ERR_FAIL_COND_MSG(!lifetime_results, "Could not init arrays.")
    memset(lifetime_results, 0, sp.particle_count * sizeof(double));

    alive_data = (bool*)malloc(sp.particle_count * sizeof(bool));
    ERR_FAIL_COND_MSG(!alive_data, "Could not init arrays.")
    memset(alive_data, false, sp.particle_count * sizeof(bool));

#ifdef RUN_CPU_SIM
    std::cout << "Simulating elastic scattering on the CPU..." << std::endl;

    QueryPerformanceCounter(&beginClock);
    if (sp.angular_speed == 0) CPUElasticScattering(sp, imp_data, lifetime_results);
    else                       CPUElasticScattering2(sp, imp_data, lifetime_results);
    QueryPerformanceCounter(&endClock);
    
    double total = 0;
    for (int i = 0; i < sp.particle_count; i++) {
        total += lifetime_results[i];
    }
    std::cout << "Average lifetime:     " << total/ sp.particle_count << " s" << std::endl;
    
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;

    for (int i = 0; i < min(sp.particle_count - 1, 100); i++)
        std::cout << lifetime_results[i] << ", ";
    return;
    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
#endif

    // Setup
    InitializeOpenCL(&ocl.deviceID, &ocl.context, &ocl.queue);

    PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    QueryPerformanceCounter(&beginClock);
    CompileOpenCLProgram(ocl.deviceID, ocl.context, source_file, &ocl.program);
    QueryPerformanceCounter(&endClock);

    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Time to build OpenCL Program: " << total_time * 1000 << " ms" << std::endl;

    PrepareOpenCLKernels();

    QueryPerformanceCounter(&beginClock);
    GPUElasticScattering(particle_count);
    QueryPerformanceCounter(&endClock);
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;

    /*std::cout.precision(64); // std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < min(particle_count - 1, 100); i++)
        std::cout << "(" << result[i] << "), ";

        */

    //std::cout << "(" << result[particle_count - 1].x << ", " << result[particle_count - 1].y << + ")" << std::endl;
}