#include <windows.h>

#include <random>
#include <limits>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"
#include "Details.h"
#include <GL/glew.h>

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
    cl_mem image;
} OCLResources;

OCLResources ocl;

typedef struct
{
    GLuint tex;
    GLuint vbo, vao;
    GLuint shader_program;

} OpenGLResources;

OpenGLResources ogl;

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

    // Shaders
    GLint success;

    std::string source = ReadShaderFile("shader.vs");
    char* vsource = new char[source.length() + 1];
    std::copy_n(source.c_str(), source.length() + 1, vsource);

    source = ReadShaderFile("shader.fs");
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

    GLuint vbo, vao;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    auto stride = 5 * sizeof(float);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    //glUseProgram(shader_program);
//    glUniform1i(glGetUniformLocation(shader_program, "texture1"), 0);


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
    
    float* pixels = new float[sp.particle_count];
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, sp.particle_row_count, sp.particle_row_count, 0, GL_RGB, GL_FLOAT, pixels);

    ocl.image = clCreateFromGLTexture2D(ocl.context, CL_MEM_WRITE_ONLY, CL_GL_OBJECT_TEXTURE2D, 0, tex, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create GLInterop texture.");
    
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
    
    clStatus = clSetKernelArg(ocl.kernel, 12, sizeof(cl_mem), (void*)&ocl.image);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");
}

void GPUElasticScattering::Compute()
{
    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    QueryPerformanceCounter(&beginClock);

    cl_int clStatus;
    size_t global_work_size[2] = { (size_t)sp.particle_row_count, (size_t)sp.particle_row_count };
    size_t local_work_size[2]  = { 10, 10 };

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus != CL_SUCCESS, clStatus, "Couldn't start kernel execution.");

    clStatus = clFinish(ocl.queue);

    /*
    std::cout << "\n\Results:" << std::endl;
    for (int i = 0; i < MIN(lifetimes.size(), 0); i++)
        std::cout << lifetimes[i] << ", ";
    std::cout << "..." << std::endl;
    */
}

void GPUElasticScattering::Draw()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);

    glUseProgram(ogl.shader_program);
    glBindVertexArray(ogl.vao);
    glDrawArrays(GL_QUADS, 0, 4);
}

GPUElasticScattering::~GPUElasticScattering()
{
    clReleaseMemObject(ocl.impb);
    clReleaseMemObject(ocl.lifetimes);
    clReleaseKernel(ocl.kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);

    glDeleteVertexArrays(1, &ogl.vao);
    glDeleteBuffers(1, &ogl.vbo);
    glDeleteProgram(ogl.shader_program);
}