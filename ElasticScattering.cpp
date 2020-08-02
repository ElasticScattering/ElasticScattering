#include <windows.h>

#include <random>
#include <limits>
#include <vector>
#include <unordered_map>
#include <assert.h>

#include "ElasticScattering.h"

typedef unsigned int uint;

const double PI  = 3.141592653589793238463;
const double PI2 = PI * 2.0;
#define RUN_CPU_SIM

inline double smod(double a, double b)
{
    return a - b * floor(a / b);
}

double ElasticScattering::GetBoundTime(const double phi, const double w, const double alpha, const bool is_electron, const bool is_future) const
{
    double remaining = smod(phi + alpha, PI * 0.5);

    double dphi;
    if (!is_electron && is_future) dphi = remaining;
    else if (!is_electron)         dphi = 2 * alpha - remaining;
    else if (is_future)            dphi = 2 * alpha - remaining;
    else                           dphi = remaining;

    return dphi / w;
}

cl_double2 ElasticScattering::GetCyclotronOrbit(const cl_double2 p, const cl_double2 velocity, const double radius, const double vf, const bool is_electron) const
{
    cl_double2 shift = { radius * velocity.y / vf, -radius * velocity.x / vf }; 

    cl_double2 center;
    if (is_electron) center = { p.x - shift.x, p.y - shift.y };
    else             center = { p.x + shift.x, p.y + shift.y };

    return center;
}

bool ElasticScattering::CirclesCross(const cl_double2 p1, const double r1, const cl_double2 p2, const double r2) const
{
    double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    if (dist_squared >= pow(r1 + r2, 2)) return false;
    if (dist_squared <= pow(r1 - r2, 2)) return false;

    return true;
}

cl_double4 ElasticScattering::GetCrossPoints(const cl_double2 p1, const double r1, const cl_double2 p2, const double r2) const
{
    const double dist_squared = pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2);
    const double dist = sqrt(dist_squared);
    const double xs = (dist_squared + r1 * r1 - r2 * r2) / (2.0 * dist);
    const double ys = sqrt(r1 * r1 - xs * xs);

    cl_double2 u = { (p2.x - p1.x) / dist, (p2.y - p1.y) / dist };

    cl_double4 z = { 
        p1.x + u.x * xs +  u.y * ys,  
        p1.y + u.y * xs + -u.x * ys, 

        p1.x + u.x * xs + -u.y * ys,
        p1.y + u.y * xs +  u.x * ys 
    };

    return z;
}

double ElasticScattering::GetPhi(const cl_double2 pos, const cl_double2 center, const double radius) const
{
    double p = (pos.x - center.x) / radius;
    assert(abs(p) < 1.0001);
    p = max(min(p, 1), -1);
    double phi = acos(p);
    
    if (pos.y < center.y)
        phi = PI2 - phi;

    return phi;
}

double ElasticScattering::GetCrossAngle(const double p, const double q, const bool clockwise) const
{
    double g = clockwise ? (p - q) : (q - p);
    return smod(g, PI2);
}

double ElasticScattering::GetCrossTime(const cl_double2 center, const cl_double2 pos, const cl_double2 ip, const double r, const double ir, const double w, const double clockwise) const
{
    const auto cross_points = GetCrossPoints(center, r, ip, ir);

    const double phi0 = GetPhi(pos, center, r);
    const double phi1 = GetPhi(cross_points.lo, center, r);
    const double phi2 = GetPhi(cross_points.hi, center, r);

    const double t1 = GetCrossAngle(phi0, phi1, clockwise) / w;
    const double t2 = GetCrossAngle(phi0, phi2, clockwise) / w;
    return min(t1, t2);
}

double ElasticScattering::CPUElasticScattering2(const cl_double2 pos, const cl_double2 vel, const SimulationParameters sp, const std::vector<cl_double2> impurities)
{
    const bool clockwise = true;

    double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
    double radius = vf / sp.angular_speed;
    auto center = GetCyclotronOrbit(pos, vel, radius, vf, clockwise);

    double lifetime = sp.particle_max_lifetime;
    for (int k = 0; k < sp.impurity_count; k++)
    {
        const cl_double2 ip = impurities[k];

        if (sp.impurity_radius_sq > pow(pos.x - ip.x, 2) + pow(pos.y - ip.y, 2))
        {
            lifetime = 0;
            break;
        }

        if (CirclesCross(center, radius, ip, sp.impurity_radius))
        {
            double t = GetCrossTime(center, pos, ip, radius, sp.impurity_radius, sp.angular_speed, clockwise);

            assert(t >= 0);
                    
            if (t < lifetime) 
                lifetime = t;
        }
    }

    return lifetime;
}

double ElasticScattering::CPUElasticScattering(const cl_double2 pos, const cl_double2 vel, const SimulationParameters sp, const std::vector<cl_double2> impurities)
{
    const cl_double2 unit = { cos(sp.phi), sin(sp.phi) };
    double lifetime = sp.particle_max_lifetime;

    for (int k = 0; k < sp.impurity_count; k++)
    {
        const cl_double2 ip = impurities[k];
        const cl_double inner = (ip.x - pos.x) * unit.x + (ip.y - pos.y) * unit.y;
        const cl_double2 projected = { pos.x + inner * unit.x, pos.y + inner * unit.y };

        const double a = pow(projected.x - ip.x, 2) + pow(projected.y - ip.y, 2);
        if (a > sp.impurity_radius_sq) {
            continue; //@Speedup, if distance is greater than current min continue as well?
        }

        double L = sqrt(sp.impurity_radius_sq - a);

        cl_double2 time_taken;
        if (vel.x != 0) time_taken = { - ((projected.x - L * unit.x) - pos.x) / vel.x, -((projected.x + L * unit.x) - pos.x) / vel.x };
        else            time_taken = { - ((projected.y - L * unit.y) - pos.y) / vel.y, -((projected.y + L * unit.y) - pos.y) / vel.y };

        if ((time_taken.s0 * time_taken.s1) < 0)
            return 0;

        if (time_taken.s0 > 0 && time_taken.s0 < lifetime) {
            lifetime = time_taken.s0;
        }
        if (time_taken.s1 > 0 && time_taken.s1 < lifetime) {
            lifetime = time_taken.s1;
        }
    }

    return lifetime;
}

void ElasticScattering::GPUElasticScattering(size_t size)
{
    cl_int clStatus;
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

void ElasticScattering::PrepareOpenCLKernels(std::vector<cl_double2> impurities, int particle_count)
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
    sp.region_size           = 1e-6;
    sp.particle_count        = 10'000; //100'000'000;
    sp.particle_row_count    = sqrt(sp.particle_count);
    sp.particle_speed        = 7e5;
    sp.particle_mass         = 5 * 9.109e-31;
    sp.impurity_count        = 100;
    sp.impurity_radius       = 1.5e-8;
    sp.impurity_radius_sq    = sp.impurity_radius * sp.impurity_radius;
    sp.alpha                 = PI / 4.0;
    sp.phi                   = sp.alpha - 1e-10;
    sp.magnetic_field        = 0;
    sp.angular_speed         = 1.602e-19 * sp.magnetic_field / sp.particle_mass;
    sp.tau                   = 1e-12;
    sp.particle_max_lifetime = sp.angular_speed == 0 ? sp.tau : min(sp.tau, GetBoundTime(sp.phi, sp.angular_speed, sp.alpha, true, false));

    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "Start region size: (" << sp.region_size << ", " << sp.region_size << ")" << std::endl;
    std::cout << "Particles:         " << sp.particle_count << std::endl;
    std::cout << "Particle speed:    " << sp.particle_speed << std::endl;
    std::cout << "Particle mass:     " << sp.particle_mass << std::endl;
    std::cout << "Impurities:        " << sp.impurity_count << std::endl;
    std::cout << "Impurity radius:   " << sp.impurity_radius << std::endl;
    std::cout << "Alpha:             " << sp.alpha << std::endl;
    std::cout << "Phi:               " << sp.phi << std::endl;
    std::cout << "Magnetic field:    " << sp.magnetic_field << std::endl;
    std::cout << "Angular speed:     " << sp.angular_speed << std::endl;
    std::cout << "Tau:               " << sp.tau << std::endl; 
    std::cout << "Particle max lifetime: " << sp.particle_max_lifetime << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    ERR_FAIL_COND_MSG(pow(sp.particle_row_count, 2) != sp.particle_count, "Particles couldn't be placed in a square grid");
    ERR_FAIL_COND_MSG(sp.alpha > (PI / 4.0), "Alpha should not be greater than pi/4.");
    ERR_FAIL_COND_MSG(sp.alpha <= 0, "Alpha should be positive.");
    ERR_FAIL_COND_MSG(sp.angular_speed < 0, "Angular speed (w) should be positive");
    ERR_FAIL_COND_MSG(sp.magnetic_field < 0, "Magnetic field strength (B) should be positive");

    // Initialize arrays.
    std::cout << "Impurity region: " << -sp.particle_speed * sp.tau << ", " << sp.region_size + sp.particle_speed * sp.tau << std::endl;
    std::uniform_real_distribution<double> unif(-sp.particle_speed * sp.tau, sp.region_size + sp.particle_speed * sp.tau);
    std::random_device r;
    std::default_random_engine re(0);

    std::vector<cl_double2> impurities(sp.impurity_count);
    for (int i = 0; i < sp.impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    lifetimes.resize(sp.particle_count, 0);

#ifdef RUN_CPU_SIM
    std::cout << "Simulating elastic scattering on the CPU..." << std::endl;

    QueryPerformanceCounter(&beginClock);
    for (int j = 0; j < sp.particle_row_count; j++)
        for (int i = 0; i < sp.particle_row_count; i++)
        {
            cl_double2 pos;
            pos.x = sp.region_size * (double(i) / sp.particle_row_count);
            pos.y = sp.region_size * (double(j) / sp.particle_row_count);

            double lifetime = sp.particle_max_lifetime;

            cl_double2 vel = { sp.particle_speed * cos(sp.phi), sp.particle_speed * sin(sp.phi) };

            double res = (sp.angular_speed == 0) ? CPUElasticScattering(pos, vel, sp, impurities) : CPUElasticScattering2(pos, vel, sp, impurities);
            lifetimes[j * sp.particle_row_count + i] = res;
        }

    QueryPerformanceCounter(&endClock);
    
    double total = 0;
    for (int i = 0; i < sp.particle_count; i++) {
        total += lifetimes[i];
    }
    std::cout << "Average lifetime:     " << total/ sp.particle_count << " s" << std::endl;
    
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;

    std::cout << "\n\nSorted results:" << std::endl;
    for (int i = 0; i < min(sp.particle_count, 200); i++)
        std::cout << lifetimes[i] << ", ";
    std::cout << "..." << std::endl;

    MakeTexture(sp);

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

    PrepareOpenCLKernels(impurities, sp.particle_count);

    QueryPerformanceCounter(&beginClock);
    GPUElasticScattering(sp.particle_count);
    QueryPerformanceCounter(&endClock);
    total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;
}

std::vector<float>  ElasticScattering::GetPixels()
{
    return pixels;
}

void ElasticScattering::MakeTexture(const SimulationParameters sp)
{
    double itau = 1.0 / sp.particle_max_lifetime;
    pixels.clear();
    pixels.resize(sp.particle_count * 3L);
    size_t j = 0;
    for (int i = 0; i < sp.particle_count; i++)
    {
        float k = float(lifetimes[i] * itau);
        if (k == 0) {
            pixels[j] = 1.0f;
            pixels[j + 1L] = 0.0f;
            pixels[j + 2L] = 0.0f;
        }
        else {
            pixels[j] = k;
            pixels[j + 1L] = k;
            pixels[j + 2L] = k;
        }
        j += 3;
    }
}