#include "windows.h"
#include "Simulation.h"
#include "src/utils/OpenCLUtils.h"

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_command_queue queue;
} OCLResources;

typedef struct
{
    cl_program program_lifetimes;
    cl_kernel lifetimes_kernel;

    cl_mem raw_lifetimes;

    cl_mem particle_settings;
    cl_mem impurity_settings;
    cl_mem simulation_settings;

    cl_mem impurities;
    cl_mem cell_indices;
    cl_mem metrics;
} OCLSimResources;

typedef struct
{
    cl_program program_integration;

    cl_kernel apply_max_lifetime;
    cl_kernel apply_sigma_comp_kernel;
    cl_kernel apply_simpson_weights;
    cl_kernel integrate_particle_kernel;
    cl_kernel sum_kernel;

    cl_mem lifetimes;

    cl_mem sigma_lifetimes_xx;
    cl_mem sigma_lifetimes_xy;

    cl_mem lifetimes_particle;

    cl_mem sigma_xx;
    cl_mem sigma_xx_particle;
    cl_mem sigma_xy;
    cl_mem sigma_xy_particle;

    cl_mem incomplete_sum; // Nog steeds nodig?
} OCLIntegrationResources;

OCLResources ocl;
OCLSimResources ocl_scatter;
OCLIntegrationResources ocl_integration;

std::vector<Sigma> SimulationCL::ComputeSigmas(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
    Metrics metrics;
    ComputeLifetimes(magnetic_field, grid, metrics);

    std::vector<Sigma> results(temperatures.size());
    for (int i = 0; i < temperatures.size(); i++)
        results[i] = DeriveTemperature(temperatures[i]);

    return results;
}

std::vector<IterationResult> SimulationCL::ComputeSigmasWithImages(const double magnetic_field, const std::vector<double>& temperatures, const Grid& grid, SampleMetrics& sample_metrics)
{
    Metrics metrics;
    ComputeLifetimes(magnetic_field, grid, metrics);

    std::vector<IterationResult> results(temperatures.size());
    for (int i = 0; i < temperatures.size(); i++)
        results[i] = DeriveTemperatureWithImages(temperatures[i]);

    return results;
}

void SimulationCL::ComputeLifetimes(const double magnetic_field, const Grid& grid, Metrics& metrics)
{
    cl_int clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.particle_settings, CL_TRUE, 0, sizeof(ParticleSettings), (void*)&ps, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.impurity_settings, CL_TRUE, 0, sizeof(ImpuritySettings), (void*)&is, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.simulation_settings, CL_TRUE, 0, sizeof(SimulationSettings), (void*)&ss, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.metrics, CL_TRUE, 0, sizeof(Metrics), (void*)&metrics, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    {
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 0, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 1, sizeof(cl_mem), (void*)&ocl_scatter.particle_settings);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 2, sizeof(cl_mem), (void*)&ocl_scatter.impurity_settings);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 3, sizeof(cl_mem), (void*)&ocl_scatter.impurities);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 4, sizeof(cl_mem), (void*)&ocl_scatter.cell_indices);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 5, sizeof(cl_mem), (void*)&ocl_scatter.raw_lifetimes);
        clStatus = clSetKernelArg(ocl_scatter.lifetimes_kernel, 6, sizeof(cl_mem), (void*)&ocl_scatter.metrics);
    }

    // Ipv WI te maken die elk een particle behandelen, laten we elk WI één phi waarde/lifetime doen.
    // @Optimize Verifiëer dat dit ook daadwerkelijk achtereenvolgende phi's berekend in één WG.
    size_t global_work_size[3] = { (size_t)ss.positions_per_row, (size_t)ss.positions_per_row, ss.particles_per_position };
    const size_t values_in_work_group = min(ss.particles_per_position, 256);
    size_t local_work_size[3] = { 256 / values_in_work_group, 1, values_in_work_group }; 

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_scatter.lifetimes_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel.");
    clFinish(ocl.queue);

    std::vector<Metrics> metrics_holder(1);
    clEnqueueReadBuffer(ocl.queue, ocl_scatter.metrics, CL_TRUE, 0, sizeof(Metrics), metrics_holder.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
    metrics = metrics_holder[0];

    std::vector<double> lifetimes(ss.total_particles);

    clEnqueueReadBuffer(ocl.queue, ocl_scatter.raw_lifetimes, CL_TRUE, 0, sizeof(Metrics), lifetimes.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back lifetimes.");
}


Sigma SimulationCL::DeriveTemperature(const double temperature) const
{
    size_t global_work_size[3] = { (size_t)ss.positions_per_row, (size_t)ss.positions_per_row, ss.particles_per_position };
    const size_t values_in_work_group = min(ss.particles_per_position, 256);
    size_t local_work_size[3] = { 256 / values_in_work_group, 1, values_in_work_group };

    cl_int clStatus;
    double tau = GetTau(temperature);
    
    // Apply max lifetime.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 0, sizeof(cl_mem), (void*)&ocl_scatter.raw_lifetimes);
        
        double default_max_lifetime = GetDefaultMaxLifetime(tau);
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 0, sizeof(double), (void*)&default_max_lifetime);
        clStatus = clSetKernelArg(ocl_integration.apply_max_lifetime, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_max_lifetime, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start apply_max_lifetime execution.");
    }

    // Apply sigma components.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.lifetimes);
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 2, sizeof(cl_mem), (void*)&ocl_scatter.particle_settings);
        clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 3, sizeof(double), (void*)&tau);

        int mode = 0;

        {
            mode = MODE_SIGMA_XX;
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 4, sizeof(int), (void*)&mode);
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 5, sizeof(cl_mem), (void*)&ocl_integration.sigma_lifetimes_xx);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_sigma_comp_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_sigma_comp_kernel xx execution.");
        }

        {
            mode = MODE_SIGMA_XY;
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 4, sizeof(int), (void*)&mode);
            clStatus = clSetKernelArg(ocl_integration.apply_sigma_comp_kernel, 5, sizeof(cl_mem), (void*)&ocl_integration.sigma_lifetimes_xy);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_sigma_comp_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_sigma_comp_kernel xy execution.");
        }
    }

    // Apply simpson weights.
    {
        clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 1, sizeof(cl_mem), (void*)&ocl_scatter.simulation_settings);
        clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 2, sizeof(int), (void*)&ss.particles_per_quadrant);

        {
            clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 0, sizeof(cl_mem), (void*)&ocl_integration.sigma_lifetimes_xx);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_simpson_weights, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_simpson_weights execution.");
        }

        {
            clStatus = clSetKernelArg(ocl_integration.apply_simpson_weights, 0, sizeof(cl_mem), (void*)&ocl_integration.sigma_lifetimes_xy);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.apply_simpson_weights, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start apply_simpson_weights execution.");
        }
    }

    Sigma sigma;

    // Sum sigma xx/xy buffers.
    {
        // @Optimize, voer de som nog vaker uit tot een aantal elementen.
        clStatus = clSetKernelArg(ocl_integration.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl_integration.incomplete_sum);
        clStatus = clSetKernelArg(ocl_integration.sum_kernel, 2, sizeof(double) * values_in_work_group, nullptr);

        const size_t half_size = ss.total_particles / 2;
        const size_t incomplete_sum_size = half_size / values_in_work_group;
        {
            clStatus = clSetKernelArg(ocl_integration.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.sigma_lifetimes_xx);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.sum_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start sum_kernel execution.");

            std::vector<double> incomplete_sum(incomplete_sum_size);
            clEnqueueReadBuffer(ocl.queue, ocl_integration.incomplete_sum, CL_TRUE, 0, sizeof(double) * incomplete_sum_size, incomplete_sum.data(), 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Failed to read back incomplete sum XX.");

            for (int i = 0; i < incomplete_sum.size(); i++)
                sigma.xx += incomplete_sum[i];
        }

        {
            clStatus = clSetKernelArg(ocl_integration.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl_integration.sigma_lifetimes_xy);

            clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl_integration.sum_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Couldn't start sum_kernel execution.");

            std::vector<double> incomplete_sum(incomplete_sum_size);
            clEnqueueReadBuffer(ocl.queue, ocl_integration.incomplete_sum, CL_TRUE, 0, sizeof(double) * incomplete_sum_size, incomplete_sum.data(), 0, nullptr, nullptr);
            CL_FAIL_CONDITION(clStatus, "Failed to read back incomplete sum XY.");

            for (int i = 0; i < incomplete_sum.size(); i++)
                sigma.xy += incomplete_sum[i];
        }
    }

    double factor = GetSigmaIntegrandFactor(tau);
    sigma.xx *= factor;
    sigma.xy *= factor;

    return sigma;
}

IterationResult SimulationCL::DeriveTemperatureWithImages(const double temperature) const
{

    IterationResult ir;

    return ir;

    
    
    /*
    // temperature is niet genoeg, moet sp meegeven...
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(settings), (void*)&sp, 0, nullptr, nullptr);
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

void SimulationCL::UploadImpurities(const Grid& grid)
{
    cl_int clStatus;

    ocl_scatter.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(v2) * grid.GetImpurities().size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create impurities buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.impurities, CL_TRUE, 0, sizeof(v2) * grid.GetImpurities().size(), grid.GetImpurities().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't upload impurities.");

    ocl_scatter.cell_indices = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(int) * grid.GetIndex().size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create cell index buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl_scatter.cell_indices, CL_TRUE, 0, sizeof(int) * grid.GetIndex().size(), grid.GetIndex().data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't upload cell indices.");
}

void SimulationCL::PrepareKernels(const Settings& ss, const size_t items_in_workgroup)
{
}

SimulationCL::SimulationCL(int p_particles_per_row, int p_values_per_quadrant) : Simulation(p_particles_per_row, p_values_per_quadrant)
{
    InitializeOpenCL(true, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (false)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    std::cout << "\nBuilding main program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "src/sim/cl/scatter.cl", &ocl_scatter.program_lifetimes);
    std::cout << "Compilation succeeded!" << std::endl;

    std::cout << "\nBuilding integration program..." << std::endl;
    CompileOpenCLProgram(ocl.deviceID, ocl.context, "src/sim/cl/integration.cl", &ocl_integration.program_integration);
    std::cout << "Compilation succeeded!" << std::endl;

    exit(0);

    //auto clStatus = clBuildProgram(program, 1, &devices[0], "-cl-mad-enable", NULL, NULL); //build the program
    //auto clStatus = clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_NUM_DEVICES, sizeof(size_t), &nb_devices, &nbread);// Return 1 devices

    size_t* np = new size_t[1]; //Create size array

    auto clStatus = clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), np, NULL);//Load in np the size of my binary

    char** bn = new char* [1]; //Create the binary array

    for (int i = 0; i < 1; i++)  bn[i] = new char[np[i]]; // I know... it's bad... but if i use new char[np], i have a segfault... :/

    clStatus = clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARIES, sizeof(unsigned char*), bn, NULL); //Load the binary itself

    /*
    printf("%s\n", bn[0]); //Print the first binary. But here, I have some curious characters

    FILE* fp;
    fopen_s(&fp, "binar.bin", "wb");

    fwrite(bn[0], sizeof(char), np[0], fp); // Save the binary, but my file stay empty
    */

    /*
    size_t number_of_binaries;
    clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &number_of_binaries, NULL);
    size_t* binaries = new size_t[number_of_binaries];
    
    auto binary = new unsigned char *[number_of_binaries];
    clGetProgramInfo(ocl_scatter.program_lifetimes, CL_PROGRAM_BINARIES, number_of_binaries, &binary, NULL);
    */


    ocl_scatter.lifetimes_kernel      = clCreateKernel(ocl_scatter.program_lifetimes, "lifetime", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create lifetime kernel.");

    ocl_integration.apply_sigma_comp_kernel = clCreateKernel(ocl_scatter.program_lifetimes, "apply_sigma_component", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create apply_sigma_component kernel.");

    ocl_integration.apply_simpson_weights = clCreateKernel(ocl_scatter.program_lifetimes, "apply_simpson_weights", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create apply_simpson_weights kernel.");

    ocl_integration.integrate_particle_kernel = clCreateKernel(ocl_scatter.program_lifetimes, "integrate_to_particle", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create integrate_to_particle kernel.");

    ocl_integration.sum_kernel = clCreateKernel(ocl_scatter.program_lifetimes, "sum", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create sum kernel.");

    ocl_scatter.particle_settings   = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ParticleSettings), nullptr, &clStatus);
    ocl_scatter.impurity_settings   = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(ImpuritySettings), nullptr, &clStatus);
    ocl_scatter.simulation_settings = clCreateBuffer(ocl.context, CL_MEM_READ_ONLY, sizeof(SimulationSettings), nullptr, &clStatus);
}

SimulationCL::~SimulationCL()
{
    clReleaseMemObject(ocl_scatter.particle_settings);
    clReleaseMemObject(ocl_scatter.impurity_settings);
    clReleaseMemObject(ocl_scatter.simulation_settings);
    clReleaseMemObject(ocl_scatter.impurities);
    clReleaseKernel(ocl_scatter.lifetimes_kernel);
    clReleaseProgram(ocl_scatter.program_lifetimes);

    clReleaseMemObject(ocl_integration.sigma_xx);
    clReleaseMemObject(ocl_integration.sigma_xy);
    clReleaseMemObject(ocl_integration.incomplete_sum);
    clReleaseKernel(ocl_integration.apply_sigma_comp_kernel);
    clReleaseKernel(ocl_integration.apply_simpson_weights);
    clReleaseKernel(ocl_integration.integrate_particle_kernel);
    clReleaseKernel(ocl_integration.sum_kernel);
    clReleaseProgram(ocl_integration.program_integration);

    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);
}
