#include <windows.h>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"

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

    // Computes a value (sigma, raw lifetime) for a list of particles.
    cl_kernel main_kernel;

    // Computes the average lifetime of all particles.
    cl_kernel add_integral_weights_kernel;
    cl_kernel sum_kernel;

    // Computes a buffer of average main_buffer, averaging many phi's for each particle.
    cl_kernel phi_integrand_kernel;

    // Transforms buffer values to image floats (0..1)
    cl_kernel tex_kernel;

    //
    // Buffers
    //
    cl_mem parameters;
    cl_mem impurities;
    cl_mem main_buffer;
    cl_mem sum_output;

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

double last_result;

std::string ReadShaderFile(const char* shader_file)
{
    std::ifstream file(shader_file);
    std::stringstream sstream;
    sstream << file.rdbuf();

    std::string contents = sstream.str();
    return contents;
}

void GPUElasticScattering::Init(bool show_info)
{
    InitializeOpenCL(&ocl.deviceID, &ocl.context, &ocl.queue);

    if (show_info)
        PrintOpenCLDeviceInfo(ocl.deviceID, ocl.context);

    QueryPerformanceFrequency(&clockFrequency);

    const char* source_file = "scatter.cl";
    QueryPerformanceCounter(&beginClock);
    CompileOpenCLProgram(ocl.deviceID, ocl.context, source_file, &ocl.program);
    QueryPerformanceCounter(&endClock);

    if (show_info) {
        double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
        std::cout << "Time to build OpenCL Program: " << total_time * 1000 << " ms" << std::endl;
    }
    
    // OpenGL context
    GLint success;

    std::string source = ReadShaderFile("src/shaders/shader.vs");
    char* vsource = new char[source.length() + 1];
    std::copy_n(source.c_str(), source.length() + 1, vsource);

    source = ReadShaderFile("src/shaders/shader.fs");
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
}

bool GPUElasticScattering::PrepareCompute(Mode p_mode, const SimulationParameters* p_sp)
{
    cl_int clStatus;

    bool first_run = (last_sp == nullptr);

    bool nothing_changed = !first_run && mode == p_mode &&
        (sp->region_size == p_sp->region_size         && sp->dim == p_sp->dim &&
         sp->particle_speed == p_sp->particle_speed   && sp->particle_mass == p_sp->particle_mass &&
         sp->impurity_count == p_sp->impurity_count   && sp->impurity_radius == p_sp->impurity_radius &&
         sp->alpha == p_sp->alpha                     && sp->phi == p_sp->phi &&
         sp->magnetic_field == p_sp->magnetic_field   && sp->tau == p_sp->tau && 
         sp->integrand_steps == p_sp->integrand_steps && sp->clockwise == p_sp->clockwise &&
         sp->region_extends == p_sp->region_extends);

    if (nothing_changed) return false;

    sp = new SimulationParameters;
    sp->region_size        = p_sp->region_size;
    sp->region_extends     = p_sp->region_extends;
    sp->dim                = p_sp->dim;
    sp->particle_speed     = p_sp->particle_speed;
    sp->particle_mass      = p_sp->particle_mass;
    sp->impurity_count     = p_sp->impurity_count;
    sp->impurity_radius    = p_sp->impurity_radius;
    sp->alpha              = p_sp->alpha;
    sp->phi                = p_sp->phi;
    sp->magnetic_field     = p_sp->magnetic_field;
    sp->tau                = p_sp->tau;
    sp->integrand_steps    = p_sp->integrand_steps;
    sp->clockwise          = p_sp->clockwise;

    sp->particle_count     = sp->dim * sp->dim;
    sp->impurity_radius_sq = sp->impurity_radius * sp->impurity_radius;
    sp->angular_speed      = E * sp->magnetic_field / sp->particle_mass;

    if (false) {
        std::cout << "\n\n+---------------------------------------------------+" << std::endl;
        std::cout << "Simulation parameters:" << std::endl;
        std::cout << "Start region size: (" << sp->region_size << ", " << sp->region_size << ")" << std::endl;
        std::cout << "Particle speed:    " << sp->particle_speed << std::endl;
        std::cout << "Particle mass:     " << sp->particle_mass << std::endl;
        std::cout << "Impurity radius:   " << sp->impurity_radius << std::endl;
        std::cout << "Alpha:             " << sp->alpha << std::endl;
        std::cout << "Phi:               " << sp->phi << std::endl;
        std::cout << "Magnetic field:    " << sp->magnetic_field << std::endl;
        std::cout << "Angular speed:     " << sp->angular_speed << std::endl;
        std::cout << "Tau:               " << sp->tau << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;

        FAIL_CONDITION(pow(sp->dim, 2) != sp->particle_count, "Particles couldn't be placed in a square grid");
        FAIL_CONDITION(sp->alpha > (PI / 4.0), "Alpha should not be greater than pi/4.");
        FAIL_CONDITION(sp->alpha <= 0, "Alpha should be positive.");
        FAIL_CONDITION(sp->angular_speed < 0, "Angular speed (w) should be positive");
        FAIL_CONDITION(sp->magnetic_field < 0, "Magnetic field strength (B) should be positive");
    }


    if (first_run || (sp->impurity_count != last_sp->impurity_count || sp->region_extends != last_sp->region_extends || sp->region_size != last_sp->region_size)) {
        //std::cout << "Impurities:        " << sp->impurity_count << std::endl;
        PrepareImpurityBuffer();
    }

    if (first_run) {
        ocl.parameters = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(SimulationParameters), nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

        ocl.tex_kernel = clCreateKernel(ocl.program, "to_texture", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    }
    
    if (first_run || sp->particle_count != last_sp->particle_count) {
        //std::cout << "Particles:         " << sp->particle_count << std::endl;

        ocl.main_buffer = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * sp->particle_count, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

        ocl.sum_output = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(double) * sp->particle_count / 2, nullptr, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create summation buffer.");

        PrepareTexKernel();
    }

    if (first_run || p_mode != mode) {
        mode = p_mode;
        const char* kernel_name = (p_mode == Mode::AVG_LIFETIME) ? "lifetime" : "sigma_xx";
        ocl.main_kernel = clCreateKernel(ocl.program, kernel_name, &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    }
    
    if (first_run) {
        ocl.add_integral_weights_kernel = clCreateKernel(ocl.program, "add_integral_weights_2d", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

        ocl.sum_kernel = clCreateKernel(ocl.program, "sum", &clStatus);
        CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");
    }

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.parameters, CL_TRUE, 0, sizeof(SimulationParameters), (void*)sp, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.main_kernel, 0, sizeof(cl_mem), (void*)&ocl.parameters);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.main_kernel, 1, sizeof(cl_mem), (void*)&ocl.impurities);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.main_kernel, 2, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.add_integral_weights_kernel, 0, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 0, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 1, sizeof(cl_mem), (void*)&ocl.sum_output);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.sum_kernel, 2, sizeof(double) * min(sp->dim, 256), nullptr); //@todo, partial sum_kernel buffer should be synced with kernel invocation / device max work group items.
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 0, sizeof(cl_mem), (void*)&ocl.main_buffer);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    double scale = (mode == Mode::AVG_LIFETIME) ? sp->tau : sp->tau*3.0;
    clStatus = clSetKernelArg(ocl.tex_kernel, 1, sizeof(double), (void*)&scale);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    clStatus = clSetKernelArg(ocl.tex_kernel, 2, sizeof(cl_mem), (void*)&ocl.image);
    CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

    return true;
}

void GPUElasticScattering::PrepareImpurityBuffer()
{
    GenerateImpurities();

    cl_int clStatus;
    ocl.impurities = clCreateBuffer(ocl.context, CL_MEM_READ_WRITE, sizeof(v2) * impurities.size(), nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

    clStatus = clEnqueueWriteBuffer(ocl.queue, ocl.impurities, CL_TRUE, 0, sizeof(v2) * impurities.size(), impurities.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");
}

void GPUElasticScattering::PrepareTexKernel()
{
    {
        float* pixels = new float[4 * sp->particle_count];
        memset(pixels, 0, 4L * sp->particle_count);

        glGenTextures(1, &ogl.tex);
        glBindTexture(GL_TEXTURE_2D, ogl.tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, sp->dim, sp->dim, 0, GL_RGBA, GL_FLOAT, pixels);
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

double GPUElasticScattering::Compute(Mode p_mode, const SimulationParameters* p_sp)
{
    bool need_update = PrepareCompute(p_mode, p_sp);
    if (!need_update) return last_result;

    QueryPerformanceCounter(&beginClock);

    size_t global_work_size[2] = { (size_t)sp->dim, (size_t)sp->dim };
    size_t local_work_size[2] = { 16, 16 };

    cl_int clStatus;

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.main_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start main kernel execution.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.tex_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start tex_kernel kernel execution.");

    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.add_integral_weights_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start add_integral_weights kernel execution.");

    const size_t half_size = sp->particle_count / 2;
    const size_t max_work_items = min(sp->dim, 256);
    clStatus = clEnqueueNDRangeKernel(ocl.queue, ocl.sum_kernel, 1, nullptr, &half_size, &max_work_items, 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Couldn't start sum_lifetimes kernel execution.");

    clStatus = clFinish(ocl.queue);

    // @Speedup, geen copy doen met een map https://downloads.ti.com/mctools/esd/docs/opencl/memory/access-model.html
    std::vector<double> results;
    results.resize(half_size / max_work_items);
    clEnqueueReadBuffer(ocl.queue, ocl.sum_output, CL_TRUE, 0, sizeof(double) * half_size / max_work_items, results.data(), 0, nullptr, nullptr);
    CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

    double result = ComputeResult(results);

    QueryPerformanceCounter(&endClock);
    double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    //std::cout << "Simulation time: " << total_time * 1000 << " ms" << std::endl;

    double kf = sp->particle_mass * sp->particle_speed / HBAR;
    double n  = kf * kf / (PI2 * C1);
    double formula = n * E * E * sp->tau / sp->particle_mass;
    
    std::cout << "\nFormula:" << formula << std::endl;
    std::cout << "Result :" << result << std::endl;
    last_sp = sp;
    last_result = result;
    return result;
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
    clReleaseMemObject(ocl.parameters); 
    clReleaseMemObject(ocl.impurities);
    clReleaseMemObject(ocl.main_buffer);
    clReleaseMemObject(ocl.image);
    clReleaseMemObject(ocl.sum_output);
    clReleaseKernel(ocl.main_kernel);
    clReleaseKernel(ocl.tex_kernel);
    clReleaseProgram(ocl.program);
    clReleaseCommandQueue(ocl.queue);
    clReleaseContext(ocl.context);

    glDeleteVertexArrays(1, &ogl.vao);
    glDeleteBuffers(1, &ogl.vbo);
    glDeleteProgram(ogl.shader_program);
}