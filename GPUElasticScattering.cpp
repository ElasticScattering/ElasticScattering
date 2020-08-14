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

    // Kernels
    cl_kernel scatter_kernel;
    cl_kernel tex_kernel;
    cl_kernel sum_lifetimes;

    // Buffers
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

int tmp_size  = 2048;
int tmp_block = 256;

std::string ReadShaderFile(const char* shader_file)
{
    std::ifstream file(shader_file);
    std::stringstream sstream;
    sstream << file.rdbuf();

    std::string contents = sstream.str();
    return contents;
}

void GPUElasticScattering::Init(SimulationParameters p_sp)
{
    sp = p_sp;

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
    PrepareScatterKernel();
    PrepareTexKernel();
    PrepareLifetimeSumKernel();
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

    if (sp.angular_speed == 0)
    {
        clStatus = clSetKernelArg(ocl.scatter_kernel, 0, sizeof(double), (void*)&sp.region_size);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 1, sizeof(double), (void*)&sp.particle_speed);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 2, sizeof(double), (void*)&sp.tau);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 3, sizeof(double), (void*)&sp.phi);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 4, sizeof(int), (void*)&sp.impurity_count);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 5, sizeof(double), (void*)&sp.impurity_radius);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 6, sizeof(cl_mem), (void*)&ocl.impurities);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 7, sizeof(cl_mem), (void*)&ocl.image);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 7, sizeof(cl_mem), (void*)&ocl.lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
    }
    else
    {
        clStatus = clSetKernelArg(ocl.scatter_kernel, 0, sizeof(double), (void*)&sp.region_size);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 1, sizeof(double), (void*)&sp.particle_speed);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 2, sizeof(double), (void*)&sp.impurity_radius);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 3, sizeof(double), (void*)&sp.tau);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 4, sizeof(double), (void*)&sp.alpha);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 5, sizeof(double), (void*)&sp.phi);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 6, sizeof(double), (void*)&sp.angular_speed);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 7, sizeof(int), (void*)&sp.impurity_count);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 8, sizeof(cl_mem), (void*)&ocl.impurities);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

        clStatus = clSetKernelArg(ocl.scatter_kernel, 9, sizeof(cl_mem), (void*)&ocl.lifetimes);
        CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
    }
}

void GPUElasticScattering::PrepareTexKernel()
{
    {
        float* pixels = new float[4 * sp.particle_count];
        memset(pixels, 0, 4L * sp.particle_count);

        glGenTextures(1, &ogl.tex);
        glBindTexture(GL_TEXTURE_2D, ogl.tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, sp.particle_row_count, sp.particle_row_count, 0, GL_RGBA, GL_FLOAT, pixels);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

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

    clStatus = clSetKernelArg(ocl.tex_kernel, 0, sizeof(cl_mem), (void*)&ocl.lifetimes);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 1, sizeof(cl_mem), (void*)&ocl.image);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::PrepareLifetimeSumKernel()
{
    cl_int clStatus;

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

void GPUElasticScattering::Compute()
{
    sp.phi += 0.1;
    if (sp.phi > PI2)
        sp.phi = 0;
    std::cout << "Phi:               " << sp.phi << std::endl;

    cl_int clStatus;

    if (sp.angular_speed == 0) clStatus = clSetKernelArg(ocl.scatter_kernel, 3, sizeof(double), (void*)&sp.phi);
    else                       clStatus = clSetKernelArg(ocl.scatter_kernel, 5, sizeof(double), (void*)&sp.phi);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    size_t global_work_size[2] = { (size_t)sp.particle_row_count, (size_t)sp.particle_row_count };
    size_t local_work_size[2] =  { 8, 8 };

    QueryPerformanceCounter(&beginClock);

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.scatter_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start kernel execution.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.tex_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start kernel execution.");

    const size_t half_size = sp.particle_count/2;
    const size_t max_work_items = 64;
    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_lifetimes, 1, nullptr, &half_size, &max_work_items, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start kernel execution.");

    clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    std::vector<double> results;
    results.resize(half_size / max_work_items);
    clEnqueueReadBuffer(ocl.queue, ocl.summed_lifetimes, CL_TRUE, 0, sizeof(double) * half_size/max_work_items, results.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back result.");
    
    double result = 0;
    for (int i = 0; i < results.size(); i++)
        result += results[i];

    std::cout << "GPU result: " << result / sp.particle_count << std::endl;
    
    /*
    std::vector<double> results2;
    results2.resize(sp.particle_count, 1);
    clEnqueueReadBuffer(ocl.queue, ocl.lifetimes, CL_TRUE, 0, sizeof(double) * sp.particle_count, results2.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

    double result2 = 0;
    for (int i = 0; i < results2.size(); i++)
        result2 += results2[i];

    std::cout << "GPU result2: " << result2 / sp.particle_count << std::endl;
    */

    QueryPerformanceCounter(&endClock);
    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;
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