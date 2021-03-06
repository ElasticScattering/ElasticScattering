/* -------------------------------------------------------------------------
    This code is part of ElasticScattering.

    Copyright(C) 2022 Elastic Scattering developers

    This program is free software : you can redistribute it and /or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.If not, see < http://www.gnu.org/licenses/>.
   ------------------------------------------------------------------------ */

#ifndef OPENCL_UTIL_H
#define OPENCL_UTIL_H

#define CL_TARGET_OPENCL_VERSION 220

#include "ErrorMacros.h"

#include <CL/cl.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

static std::string CLErrorString(int err) {
    switch (err) {
        case CL_SUCCESS:                            return "Success!";
        case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
        case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
        case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
        case CL_OUT_OF_RESOURCES:                   return "Out of resources";
        case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
        case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
        case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
        case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
        case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
        case CL_MAP_FAILURE:                        return "Map failure";
        case CL_INVALID_VALUE:                      return "Invalid value";
        case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
        case CL_INVALID_PLATFORM:                   return "Invalid platform";
        case CL_INVALID_DEVICE:                     return "Invalid device";
        case CL_INVALID_CONTEXT:                    return "Invalid context";
        case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
        case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
        case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
        case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
        case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
        case CL_INVALID_SAMPLER:                    return "Invalid sampler";
        case CL_INVALID_BINARY:                     return "Invalid binary";
        case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
        case CL_INVALID_PROGRAM:                    return "Invalid program";
        case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
        case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
        case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
        case CL_INVALID_KERNEL:                     return "Invalid kernel";
        case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
        case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
        case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
        case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
        case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
        case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
        case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
        case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
        case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
        case CL_INVALID_EVENT:                      return "Invalid event";
        case CL_INVALID_OPERATION:                  return "Invalid operation";
        case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
        case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
        case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
        default: return "Unknown error";
    }
}

static void InitializeOpenCL(bool prefer_gpu, cl_device_id *p_deviceID, cl_context *p_ctx, cl_command_queue *p_queue)
{
    *p_deviceID = nullptr;
    *p_ctx      = nullptr;

    cl_uint num_platforms = 0;
    cl_int  cl_status = clGetPlatformIDs(0, nullptr, &num_platforms);
    CL_FAIL_CONDITION(cl_status, "No platforms found.");
    FAIL_CONDITION(num_platforms == 0, "No platforms found.");

    char             platform_vendor[256];
    cl_context       ctx = nullptr;
    cl_command_queue queue = nullptr;

    std::vector<cl_platform_id> platform_ids(num_platforms);
    clGetPlatformIDs(num_platforms, platform_ids.data(), nullptr);

    float          best_device_score = 0;
    cl_platform_id selected_platform = nullptr;
    cl_device_id   selected_device = nullptr;
     
    for (int i = 0; i < num_platforms; i++)
    {
        float score = 1.0f;
        cl_platform_id platform_id = platform_ids[i];

        cl_status = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, sizeof(platform_vendor), platform_vendor, nullptr);
        CL_FAIL_CONDITION(cl_status, "Could not get platform info.");

        bool is_intel = strncmp("Intel", platform_vendor, strlen("Intel")) == 0;
        if (prefer_gpu) 
            score *= is_intel ? 0.1f : 1;
        else 
            score *= is_intel ? 1 : 0.1f;

        cl_device_id device_id = nullptr;
        cl_device_type device_type = CL_DEVICE_TYPE_GPU;
        cl_status = clGetDeviceIDs(platform_id, device_type, 1, &device_id, nullptr);
        if (cl_status != CL_SUCCESS) {
            device_type = CL_DEVICE_TYPE_CPU;
            cl_status = clGetDeviceIDs(platform_id, device_type, 1, &device_id, nullptr);
        }

        score *= device_type == CL_DEVICE_TYPE_GPU ? 1 : 0.6f;

        if (score > best_device_score) {
            selected_device = device_id;
            best_device_score = score;
            selected_platform = platform_id;
        }
    }

    cl_context_properties props[] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)selected_platform,
        0
    };

    ctx = clCreateContext(props, 1, &selected_device, nullptr, nullptr, &cl_status);
    CL_FAIL_CONDITION(cl_status, "Couldn't create context.");

    queue = clCreateCommandQueueWithProperties(ctx, selected_device, 0, &cl_status);
    CL_FAIL_CONDITION(cl_status, "Couldn't create command queue.");

    *p_deviceID = selected_device;
    *p_ctx = ctx;
    *p_queue = queue;
}

static void CompileOpenCLProgram(const cl_device_id p_device_id, const cl_context p_ocl_context, const char *p_soure_file, cl_program *p_ocl_program)
{
    cl_int clStatus;
    
    std::string path(p_soure_file);
    std::fstream kernelFile(path);
    std::string content((std::istreambuf_iterator<char>(kernelFile)), std::istreambuf_iterator<char>());
    const char *code = new char[content.size()];
    code = content.c_str();

    cl_program program = clCreateProgramWithSource(p_ocl_context, 1, (const char**)&code, nullptr, &clStatus);
    CL_FAIL_CONDITION(clStatus, "Couldn't create program.");

    clStatus = clBuildProgram(program, 1, &p_device_id, "-D DEVICE_PROGRAM", NULL, NULL);

    if (clStatus != CL_SUCCESS) {
        size_t len;
        char buffer[4086];

        clGetProgramBuildInfo(program, p_device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        CL_FAIL_CONDITION(clStatus, buffer);
    }

    *p_ocl_program = program;
}

static void PrintOpenCLDeviceInfo(const cl_device_id device_id, const cl_context contextHdl)
{
    char        device_vendor[512];
    char        device_name[512];
    char        driver_version[512];
    char        device_profile[512];
    char        device_version[512];
    char        oclc_version[512];
    cl_uint     max_compute_units = 0;
    cl_uint     max_work_items_dim = 0;
    size_t      max_work_item_sizes[3];
    size_t      max_work_group_size;
    cl_uint     max_device_frequency;
    cl_uint     max_num_samplers = 0;
    cl_uint     min_base_addr_align_size_bits = 0;
    cl_ulong    local_mem_size;
    cl_ulong    global_mem_size;
    cl_long	    max_mem_alloc_size;
    cl_long	    max_constant_buffer_size;

    size_t      num_bytes = 0;
    cl_int      cl_status;

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, sizeof(device_vendor), &device_vendor, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(device_name), &device_name, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DRIVER_VERSION, sizeof(driver_version), &driver_version, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_PROFILE, sizeof(device_profile), &device_profile, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_VERSION, sizeof(device_version), &device_version, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_OPENCL_C_VERSION, sizeof(oclc_version), &oclc_version, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &max_work_items_dim, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(max_work_item_sizes), &max_work_item_sizes, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(cl_uint), &min_base_addr_align_size_bits, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &max_device_frequency, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_mem_size, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_long), &max_mem_alloc_size, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(cl_long), &max_constant_buffer_size, &num_bytes);
    CL_FAIL_CONDITION(cl_status, "clGetDeviceInfo() query failed.");

    size_t info_size;
    auto p_info = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, 0, nullptr, &info_size);
    char* info = (char*)_malloca(sizeof(char) * info_size);
    clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, info_size, info, nullptr);

    std::string infos(info);
    std::stringstream ss(infos);
    std::string ext;
    std::vector<std::string> extensions;

    while (std::getline(ss, ext, ' ')) {
        extensions.push_back(ext);
    }

    bool glinterop  = std::find(extensions.begin(), extensions.end(), "cl_khr_gl_sharing") != extensions.end();
    bool double_sup = std::find(extensions.begin(), extensions.end(), "cl_khr_fp64") != extensions.end();


    std::cout << "OpenCL device info:"        << std::endl;
    std::cout << "Device vendor:            " << device_vendor << std::endl;
    std::cout << "Device name:              " << device_name << std::endl;
    std::cout << "Driver version:           " << driver_version << std::endl;
    std::cout << "Device profile:           " << device_profile << std::endl;
    std::cout << "Device version:           " << device_version << std::endl;
    std::cout << "OpenCL C version:         " << oclc_version << std::endl;
    std::cout << std::endl;
    std::cout << "Max compute units:        " << max_compute_units << std::endl;
    std::cout << "Max work item dimension:  " << max_work_items_dim << std::endl;
    std::cout << "Max work item sizes:      " << max_work_item_sizes[0] << ", " << max_work_item_sizes[1] << ", " << max_work_item_sizes[2] << std::endl;
    std::cout << "Max work group size:      " << max_work_group_size << std::endl;
    std::cout << "Largest builtin type:     " << min_base_addr_align_size_bits << " (bits)" << std::endl;
    std::cout << "Local memory size:        " << (float)local_mem_size << std::endl;
    std::cout << "Global memory size:       " << (float)global_mem_size << std::endl;
    std::cout << "Max memory alloc size:    " << (float)max_mem_alloc_size << std::endl;
    std::cout << "Max constant buffer size: " << max_constant_buffer_size << std::endl;
    std::cout << "Max clock frequency:      " << max_device_frequency << std::endl;
    std::cout << "Double Precision:         " << (double_sup ? "True" : "False") << std::endl;
    //std::cout << "GL Interop:                     " << (glinterop ? "True" : "False") << std::endl;
    //std::cout << "Device image2d max width:       " << uMaxImage2DWidth << std::endl;
}

#endif // OPENCL_UTIL_H