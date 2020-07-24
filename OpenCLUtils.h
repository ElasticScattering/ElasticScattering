#ifndef OPENCL_UTIL_H
#define OPENCL_UTIL_H

#include <CL/cl.h>
#include <CL/cl_gl.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "glew.h"
#include <wglew.h>
#include <cl/cl_gl_ext.h>

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

#define CL_ERR_FAIL_COND_MSG(err, msg)                                                                                                                                                          \
    if (err != CL_SUCCESS) {                                                                                                                                                                    \
        std::cout << "\x1B[31mError [" << CLErrorString(err) << "]\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl;   \
        system("pause");                                                                                                                                          \
        exit(0);                                                                                                                                                                                \
    }

#define ERR_FAIL_COND_MSG(cond, msg)                                                                                                                            \
	if (cond) {                                                                                                                                                 \
		std::cout << "\x1B[31mError\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl;  \
        system("pause");                                                                                                                                          \
		exit(0);                                                                                                                                                \
	}

static void InitializeOpenCL(cl_device_id *p_deviceID, cl_context *p_ctx, cl_command_queue *p_queue) 
{
    *p_deviceID = nullptr;
    *p_ctx      = nullptr;

    cl_uint num_platforms = 0;
    cl_int  cl_status = clGetPlatformIDs(0, nullptr, &num_platforms);
    CL_ERR_FAIL_COND_MSG(cl_status, "No platforms found.");
    ERR_FAIL_COND_MSG(num_platforms == 0, "No platforms found.");

    char             platform_vendor[256];
    char             device_version[256];
    char             lang_version[256];
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
        CL_ERR_FAIL_COND_MSG(cl_status, "Could not get platform info.");

        bool is_intel = strncmp("Intel", platform_vendor, strlen("Intel")) == 0;
        score *= is_intel ? 0.1f : 1;

        cl_device_id device_id = nullptr;
        cl_device_type device_type = CL_DEVICE_TYPE_GPU;
        cl_status = clGetDeviceIDs(platform_id, device_type, 1, &device_id, nullptr);
        if (cl_status != CL_SUCCESS) {
            device_type = CL_DEVICE_TYPE_CPU;
            cl_status = clGetDeviceIDs(platform_id, device_type, 1, &device_id, nullptr);
        }

        score *= device_type == CL_DEVICE_TYPE_GPU ? 1 : 0.3f;

        if (score > best_device_score) {
            selected_device = device_id;
            best_device_score = score;
            selected_platform = platform_id;
        }
    }

    /*
    cl_context_properties props[] =
    {
        CL_GL_CONTEXT_KHR,
        (cl_context_properties)wglGetCurrentContext(),
        CL_WGL_HDC_KHR,
        (cl_context_properties)wglGetCurrentDC(),
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)selected_platform,
        0
    };
    */

    ctx = clCreateContext(0, 1, &selected_device, nullptr, nullptr, &cl_status);
    CL_ERR_FAIL_COND_MSG(!ctx, cl_status, "Couldn't create context.");

    queue = clCreateCommandQueueWithProperties(ctx, selected_device, 0, &cl_status);
    CL_ERR_FAIL_COND_MSG(!queue, cl_status, "Couldn't create command queue.");

    *p_deviceID = selected_device;
    *p_ctx = ctx;
    *p_queue = queue;
}

static void CompileOpenCLProgram(const cl_device_id p_device_id, const cl_context p_ocl_context, const char * p_soure_file, cl_program *p_ocl_program)
{
    cl_int clStatus;
    
    std::fstream kernelFile(p_soure_file);
    std::string content((std::istreambuf_iterator<char>(kernelFile)), std::istreambuf_iterator<char>());
    const char *code = new char[content.size()];
    code = content.c_str();

    cl_program program = clCreateProgramWithSource(p_ocl_context, 1, (const char**)&code, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(clStatus, "Couldn't create program.");

    clStatus = clBuildProgram(program, 1, &p_device_id, NULL, NULL, NULL);
    if (clStatus != CL_SUCCESS) {
        size_t len;
        char buffer[2048];

        clGetProgramBuildInfo(program, p_device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        CL_ERR_FAIL_COND_MSG(clStatus, buffer);
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
    cl_uint     min_base_addr_align_size_bytes = 0;
    cl_ulong    local_mem_size;
    cl_long	    max_mem_alloc_size;

    size_t      num_bytes = 0;
    cl_int      cl_status;

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, sizeof(device_vendor), &device_vendor, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(device_name), &device_name, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DRIVER_VERSION, sizeof(driver_version), &driver_version, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_PROFILE, sizeof(device_profile), &device_profile, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_VERSION, sizeof(device_version), &device_version, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_OPENCL_C_VERSION, sizeof(oclc_version), &oclc_version, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &max_work_items_dim, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(max_work_item_sizes), &max_work_item_sizes, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(cl_uint), &min_base_addr_align_size_bits, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, sizeof(cl_uint), &min_base_addr_align_size_bytes, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");
    
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &max_device_frequency, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_long), &max_mem_alloc_size, &num_bytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");

    size_t info_size;
    auto p_info = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, 0, nullptr, &info_size);
    char* info = (char*)_malloca(sizeof(char) * info_size);
    clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, info_size, info, nullptr);

    /* Doesn't work
    cl_uint	uMaxImage2DWidth;
    cl_status = clGetDeviceInfo(device_id, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(cl_uint), &uMaxImage2DWidth, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(cl_status, "clGetDeviceInfo() query failed.");
    */

    std::cout << "OpenCL device info:"              << std::endl;
    std::cout << "Device vendor:                  " << device_vendor << std::endl;
    std::cout << "Device name:                    " << device_name << std::endl;
    std::cout << "Driver version:                 " << driver_version << std::endl;
    std::cout << "Device profile:                 " << device_profile << std::endl;
    std::cout << "Device version:                 " << device_version << std::endl;
    std::cout << "Device OpenCL C version:        " << oclc_version << std::endl;
    std::cout << "Device max compute units:       " << max_compute_units << std::endl;
    std::cout << "Device max work item dim.:      " << max_work_items_dim << std::endl;
    std::cout << "Device max work item sizes:     " << max_work_item_sizes[0] << ", " << max_work_item_sizes[1] << ", " << max_work_item_sizes[2] << std::endl;
    std::cout << "Device max work group size:     " << max_work_group_size << std::endl;
    std::cout << "Device mem. base addr. align:   " << min_base_addr_align_size_bits << " (bits)" << std::endl;
    std::cout << "Device min datatype align size: " << min_base_addr_align_size_bytes << " (bytes)" << std::endl;
    std::cout << "Device max clock frequency:     " << max_device_frequency << std::endl;
    std::cout << "Device local mem. size:         " << (float)local_mem_size << std::endl;
    std::cout << "Device max mem alloc size:      " << (float)max_mem_alloc_size << std::endl;
    std::cout << "Device extensions:              " << info << std::endl;
    //std::cout << "Device image2d max width:       " << uMaxImage2DWidth << std::endl;
}

#endif // OPENCL_UTIL_H