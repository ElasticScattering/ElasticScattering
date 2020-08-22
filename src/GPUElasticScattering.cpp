#include <windows.h>

#include <random>
#include <limits>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"
#include "Details.h"
#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>

typedef struct
{
    cl_device_id deviceID;
    cl_context context;
    cl_program program;
    cl_command_queue queue;

    //
    // Kernels
    //

    // Computes a buffer of raw lifetimes.
    cl_kernel scatter_kernel;

    // Computes the average lifetime of all particles.
    cl_kernel add_integral_weights;
    cl_kernel sum_lifetimes;

    // Computes a buffer of average lifetimes, averaging many phi's for each particle.
    cl_kernel phi_integrand;

    // Transforms buffer values to image floats (0..1)
    cl_kernel tex_kernel;

    //
    // Buffers
    //
    cl_mem impurities;
    cl_mem lifetimes;
    cl_mem summed_lifetimes;
    cl_mem sigma_xx;

    cl_mem image;
} OCLResources;

typedef struct
{
    GLuint tex, tex2;
    GLuint vbo, vao;
    GLuint shader_program;
} OpenGLResources;

OCLResources ocl;
OpenGLResources ogl;

LARGE_INTEGER beginClock, endClock, clockFrequency;

std::string ReadShaderFile(const char* shader_file)
{
    std::ifstream file(shader_file);
    std::stringstream sstream;
    sstream << file.rdbuf();

    std::string contents = sstream.str();
    return contents;
}

void GPUElasticScattering::Init(Mode p_mode, SimulationParameters p_sp)
{
    sp = p_sp;
    mode = p_mode;

    QueryPerformanceFrequency(&clockFrequency);
    
    char		deviceStr[256];
    char		vendorStr[256];
    
    InitializeOpenCL(&ocl.deviceID, &ocl.context, &ocl.queue);
    PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    const char* source_file = (sp.angular_speed == 0) ? "scatter0.cl" : "scatterB.cl";
    QueryPerformanceCounter(&beginClock);
    CompileOpenCLProgram(ocl.deviceID, ocl.context, source_file, &ocl.program);
    QueryPerformanceCounter(&endClock);

    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Time to build OpenCL Program: " << total_time * 1000 << " ms" << std::endl;

    // OpenGL context
    GLint success;

    std::string source = ReadShaderFile("shaders/shader.vs");
    char* vsource = new char[source.length() + 1];
    std::copy_n(source.c_str(), source.length() + 1, vsource);

    source = ReadShaderFile("shaders/shader.fs");
    char* fsource = new char[source.length() + 1];
    std::copy_n(source.c_str(), source.length() + 1, fsource);

    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, &vsource, nullptr);
    glCompileShader(vert_shader);

    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[2048];
        glGetShaderInfoLog(vert_shader, 2048, nullptr, infoLog);
        std::cout << "Failed to compile vertex shader. Info:\n" << infoLog << std::endl;
    }

    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(frag_shader, 1, &fsource, nullptr);
    glCompileShader(frag_shader);

    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[2048];
        glGetShaderInfoLog(frag_shader, 2048, nullptr, infoLog);
        std::cout << "Failed to compile fragment shader. Info:\n" << infoLog << std::endl;
    }

    ogl.shader_program = glCreateProgram();
    glAttachShader(ogl.shader_program, vert_shader);
    glAttachShader(ogl.shader_program, frag_shader);
    glLinkProgram(ogl.shader_program);
    glUseProgram(ogl.shader_program);
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader);

    // Texture
    float vertices[] =
    {
        // Position,        Tex coord
        -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
         1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
         1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
        -1.0f,  1.0f, 0.0f, 0.0f, 1.0f
    };

    glGenVertexArrays(1, &ogl.vao);
    glGenBuffers(1, &ogl.vbo);

    glBindVertexArray(ogl.vao);

    glBindBuffer(GL_ARRAY_BUFFER, ogl.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    auto stride = 5 * sizeof(float);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    impurities.clear();
    impurities.resize(sp.impurity_count);

    std::cout << "Impurity region: " << -sp.particle_speed * sp.tau << ", " << sp.region_size + sp.particle_speed * sp.tau << std::endl;
    std::uniform_real_distribution<double> unif(-sp.particle_speed * sp.tau, sp.region_size + sp.particle_speed * sp.tau);
    std::random_device r;
    std::default_random_engine re(0);

    for (int i = 0; i < sp.impurity_count; i++)
        impurities[i] = { unif(re), unif(re) };

    PrepareOpenCLKernels();
}

void GPUElasticScattering::PrepareOpenCLKernels()
{
    if (mode == Mode::SIGMA_XX)
    {
        PrepareIntegrandKernel();
    }
    else if (mode == Mode::AVG_LIFETIME)
    {
        PrepareScatterKernel();
        PrepareLifetimeSumKernel();
    }
    PrepareTexKernel();
}

void GPUElasticScattering::PrepareScatterKernel()
{
    cl_int clStatus;

    ocl.scatter_kernel = clCreateKernel(ocl.program, "lifetime", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(v2) * impurities.size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * impurities.size(), impurities.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

    ocl.lifetimes = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * sp.particle_count, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 0, sizeof(SimulationParameters), (void*)&sp);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 1, sizeof(cl_mem), (void*)&ocl.impurities);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.scatter_kernel, 2, sizeof(cl_mem), (void*)&ocl.lifetimes);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::PrepareTexKernel()
{
    {
        float* pixels = new float[4 * sp.particle_count];
        memset(pixels, 0, 4L * sp.particle_count);

        glGenTextures(1, &ogl.tex);
        glBindTexture(GL_TEXTURE_2D, ogl.tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, sp.dim, sp.dim, 0, GL_RGBA, GL_FLOAT, pixels);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); //GL_LINEAR
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        glUseProgram(ogl.shader_program);
        glUniform1i(glGetUniformLocation(ogl.shader_program, "texture1"), 0);
        delete[] pixels;
    }

    cl_int clStatus; 

    ocl.tex_kernel = clCreateKernel(ocl.program, "to_texture", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    
    ocl.image = clCreateFromGLTexture(ocl.context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, ogl.tex, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create GLInterop texture.");

    glUseProgram(ogl.shader_program);
    glUniform1i(glGetUniformLocation(ogl.shader_program, "texture1"), 0);

    auto buffer = (mode == Mode::SIGMA_XX) ? ocl.sigma_xx : ocl.lifetimes;

    clStatus = clSetKernelArg(ocl.tex_kernel, 0, sizeof(cl_mem), (void*)&buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 1, sizeof(double), (void*)&sp.tau);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 2, sizeof(cl_mem), (void*)&ocl.image);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::PrepareLifetimeSumKernel()
{
    cl_int clStatus;

    ocl.add_integral_weights = clCreateKernel(ocl.program, "add_integral_weights_2d", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    clStatus = clSetKernelArg(ocl.add_integral_weights, 0, sizeof(cl_mem), (void*)&ocl.lifetimes);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    //////////////
    
    ocl.sum_lifetimes = clCreateKernel(ocl.program, "sum", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.summed_lifetimes = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * sp.particle_count/2, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create summation buffer.");
    
    clStatus = clSetKernelArg(ocl.sum_lifetimes, 0, sizeof(cl_mem), (void*)&ocl.lifetimes);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_lifetimes, 1, sizeof(cl_mem), (void*)&ocl.summed_lifetimes);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_lifetimes, 2, sizeof(double) * 128, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::PrepareIntegrandKernel()
{
    cl_int clStatus;

    ocl.phi_integrand = clCreateKernel(ocl.program, "sigma_xx", &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

    ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(v2) * impurities.size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * impurities.size(), impurities.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

    ocl.sigma_xx = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * sp.particle_count, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

    clStatus = clSetKernelArg(ocl.phi_integrand, 0, sizeof(SimulationParameters), (void*)&sp);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.phi_integrand, 2, sizeof(cl_mem), (void*)&ocl.impurities);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.phi_integrand, 3, sizeof(cl_mem), (void*)&ocl.sigma_xx);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

double GPUElasticScattering::Compute()
{
    QueryPerformanceCounter(&beginClock);

    size_t global_work_size[2] = { (size_t)sp.dim, (size_t)sp.dim };
    size_t local_work_size[2]  = { 16, 16 };

    if (mode == Mode::AVG_LIFETIME)
    {
        sp.phi += 0.1;
        if (sp.phi > PI2)
            sp.phi = 0;
        std::cout << "Phi:               " << sp.phi << std::endl;

        cl_int clStatus;
        
        clStatus = clSetKernelArg(ocl.scatter_kernel, 0, sizeof(SimulationParameters), (void*)&sp);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.scatter_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start scatter_kernel kernel execution.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.tex_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start tex_kernel kernel execution.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.add_integral_weights, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

        const size_t half_size = sp.particle_count / 2;
        const size_t max_work_items = MIN(sp.dim, 256);
        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_lifetimes, 1, nullptr, &half_size, &max_work_items, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start sum_lifetimes kernel execution.");

        clStatus = clFinish(ocl.queue);

        // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
        std::vector<double> results;
        results.resize(half_size / max_work_items);
        clEnqueueReadBuffer(ocl.queue, ocl.summed_lifetimes, CL_TRUE, 0, sizeof(double) * half_size / max_work_items, results.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

        double total = 0;
        for (int i = 0; i < results.size(); i++)
        {
            std::cout << results[i] << ", ";
            total += results[i];
        }

        double simulated_particle_count = (sp.dim - 1) * (sp.dim - 1);
        double result = total / simulated_particle_count;
        std::cout << "GPU result: " << result << std::endl;
        return result;
    }
    else if (mode == Mode::SIGMA_XX) {
        cl_int clStatus;

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.phi_integrand, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start kernel execution.");

        clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.tex_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Couldn't start kernel execution.");

        clStatus = clFinish(ocl.queue);
        
        /*
        std::vector<double> results;
        results.resize(sp.particle_count);
        clEnqueueReadBuffer(ocl.queue, ocl.sigma_xx, CL_TRUE, 0, sizeof(double) * sp.particle_count, results.data(), 0, nullptr, nullptr);
        CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

        for (int i = 0; i < MIN(600, results.size()); i++)
            std::cout << results[i] << ", ";
        */
    }
   
    QueryPerformanceCounter(&endClock);
    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;

    return 0;
}

void GPUElasticScattering::Draw()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ogl.tex);

    glUseProgram(ogl.shader_program);
    glBindVertexArray(ogl.vao);
    glDrawArrays(GL_QUADS, 0, 4);
}

GPUElasticScattering::~GPUElasticScattering()
{
    clReleaseMemObject(ocl.impurities);
    clReleaseMemObject(ocl.lifetimes);
    clReleaseMemObject(ocl.image);
    clReleaseMemObject(ocl.summed_lifetimes);
    clReleaseKernel(ocl.scatter_kernel);
    clReleaseKernel(ocl.tex_kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);

    glDeleteVertexArrays(1, &ogl.vao);
    glDeleteBuffers(1, &ogl.vbo);
    glDeleteProgram(ogl.shader_program);
}