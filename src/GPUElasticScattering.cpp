#include <windows.h>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"

#ifndef NO_WINDOW
#include "utils/OpenGLUtils.h"
#endif

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;

    // Computes a value (sigma, raw lifetime) for a list of particles.
    cl_kernel scatter_kernel;

    // Computes the average lifetime of all particles.
    cl_kernel add_integral_weights_kernel;
    cl_kernel sum_kernel;

    // Transforms buffer values to image floats (0..1)
    cl_kernel tex_kernel;

    // Buffers
    cl_mem parameters;
    cl_mem impurities;
    cl_mem main_buffer;
    cl_mem sum_output;

    cl_mem image;
} OCLResources;

OCLResources ocl;

bool GPUElasticScattering::Compute(ScatteringParameters& p_sp, double &result)
{
    bool need_update = PrepareCompute(p_sp);
    if (!need_update) return false;

    size_t global_work_size[2] = { (size_t)sp.dim, (size_t)sp.dim };
    size_t local_work_size[2] = { min(sp.dim, 256), 256 / min(sp.dim, 256) };

    cl_int clStatus;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.scatter_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel execution.");

#ifndef NO_WINDOW
    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.tex_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start tex_kernel kernel execution.");
#endif //NO_WINDOW

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.add_integral_weights_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

    const size_t half_size = particle_count / 2;
    const size_t max_work_items = min(sp.dim, 256);
    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_kernel, 1, nullptr, &half_size, &max_work_items, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start sum_lifetimes kernel execution.");

    clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    std::vector<double> results;
    results.resize(half_size / max_work_items);
    clEnqueueReadBuffer(ocl.queue, ocl.sum_output, CL_TRUE, 0, sizeof(double) * half_size / max_work_items, results.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

    result = ComputeResult(results);

    return true;
}

bool GPUElasticScattering::PrepareCompute(ScatteringParameters &p_sp)
{
    cl_int clStatus;

    CompleteSimulationParameters(p_sp);

    if (!first_run && !AnythingChanged(p_sp)) return false;

    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed  = (sp.dim != p_sp.dim);
    
    sp = p_sp;

    if (first_run || impurities_changed) {
        GenerateImpurities(sp);

        ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(v2) * impurities.size(), nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

        clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * impurities.size(), impurities.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");
    }

    if (first_run) {
        ocl.parameters = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(ScatteringParameters), nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");
#ifndef NO_WINDOW
        ocl.tex_kernel = clCreateKernel(ocl.program, "to_texture", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
#endif //NO_WINDOW
    }
    
    if (first_run || work_size_changed) {
        ocl.main_buffer = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * particle_count, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

        ocl.sum_output = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * particle_count / 2, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create summation buffer.");

#ifndef NO_WINDOW
        PrepareTexKernel(sp.dim);
#endif //NO_WINDOW
    }

    if (first_run) {
        ocl.scatter_kernel = clCreateKernel(ocl.program, "lifetime", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    
        ocl.add_integral_weights_kernel = clCreateKernel(ocl.program, "add_integral_weights_2d", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

        ocl.sum_kernel = clCreateKernel(ocl.program, "sum", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    }

    first_run = false;

    //clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(cl_mem), (void*)&sp, 0, nullptr, nullptr);
    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(ScatteringParameters), (void*)&sp, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 0, sizeof(cl_mem), (void*)&ocl.parameters);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 1, sizeof(cl_mem), (void*)&ocl.impurities);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 2, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.add_integral_weights_kernel, 0, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl.sum_output);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 2, sizeof(double) * min(sp.dim, 256), nullptr); //@todo, partial sum_kernel buffer should be synced with kernel invocation / device max work group items.
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

#ifndef NO_WINDOW
    clStatus = clSetKernelArg(ocl.tex_kernel, 0, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 1, sizeof(int), (void*)&sp.mode);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    //double scale = ShouldComputeSigma(sp.mode) ? sp.tau : sp.tau*3.0;
    clStatus = clSetKernelArg(ocl.tex_kernel, 2, sizeof(double), (void*)&sp.tau);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 3, sizeof(cl_mem), (void*)&ocl.image);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
#endif //NO_WINDOW

    return true;
}

#ifndef NO_WINDOW
void GPUElasticScattering::PrepareTexKernel(int dim)
{
    {
        int pixel_count = dim * dim;
        float* pixels = new float[4 * pixel_count];
        memset(pixels, 0, 4L * pixel_count);

        glGenTextures(1, &ogl.tex);
        glBindTexture(GL_TEXTURE_2D, ogl.tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, dim, dim, 0, GL_RGBA, GL_FLOAT, pixels);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); //GL_LINEAR
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        glUseProgram(ogl.shader_program);
        glUniform1i(glGetUniformLocation(ogl.shader_program, "texture1"), 0);
        delete[] pixels;
    }

    cl_int clStatus; 

    ocl.image = clCreateFromGLTexture(ocl.context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, ogl.tex, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create GLInterop texture.");

    glUseProgram(ogl.shader_program);
    glUniform1i(glGetUniformLocation(ogl.shader_program, "texture1"), 0);
}
#endif // NO_WINDOW

GPUElasticScattering::GPUElasticScattering()
{
    InitParameters init;
    init.dont_show_info = false;
    init.use_gpu = true;

    InitializeOpenCL(init.use_gpu, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (!init.dont_show_info)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    CompileOpenCLProgram(ocl.deviceID, ocl.context, "scatter.cl", &ocl.program);

#ifndef NO_WINDOW
    OpenGLUtils o;
    o.Init(ogl.vbo, ogl.vao, ogl.shader_program);
#endif //NO_WINDOW
}

GPUElasticScattering::GPUElasticScattering(const InitParameters& init)
{
    InitializeOpenCL(init.use_gpu, &ocl.deviceID, &ocl.context, &ocl.queue);
    if (!init.dont_show_info)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    CompileOpenCLProgram(ocl.deviceID, ocl.context, "scatter.cl", &ocl.program);
    
#ifndef NO_WINDOW
    OpenGLUtils o;
    o.Init(ogl.vbo, ogl.vao, ogl.shader_program);
#endif //NO_WINDOW
}

GPUElasticScattering::~GPUElasticScattering()
{
    clReleaseMemObject(ocl.parameters); 
    clReleaseMemObject(ocl.impurities);
    clReleaseMemObject(ocl.main_buffer);
    clReleaseMemObject(ocl.image);
    clReleaseMemObject(ocl.sum_output);
    clReleaseKernel(ocl.scatter_kernel);
    clReleaseKernel(ocl.tex_kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);

#ifndef NO_WINDOW
    glDeleteVertexArrays(1, &ogl.vao);
    glDeleteBuffers(1, &ogl.vbo);
    glDeleteProgram(ogl.shader_program);
#endif
}
