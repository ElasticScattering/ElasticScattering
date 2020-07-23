#ifndef UTIL_H
#define UTIL_H

#include <CL/cl.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

std::string CLErrorString(int err) {
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

#define CL_ERR_FAIL_COND_MSG(err, msg)                                                                                                         \
    if (err != CL_SUCCESS) {                                                                                                                                          \
        std::cout << "\x1B[31mError [" << CLErrorString(err) << "]\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl; \
        exit(0);                                                                                                                                     \
    } else                                                                                                                                               \
        ((void)0);

#define ERR_FAIL_COND_MSG(cond, msg)                                                                               \
	if (cond) {                                                                                           \
		std::cout << "\x1B[31mError\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << "\n\tin function: " << __FUNCTION__ << std::endl; \
		exit(0);                                                                                        \
	} else                                                                                                \
		((void)0);

void InitializeOpenCL(char *p_device_str, char *p_vender_str, cl_device_id *p_deviceID, cl_context *p_ctx, cl_command_queue *p_queue) 
{
    *p_deviceID = nullptr;
    *p_ctx      = nullptr;

    cl_int clStatus;
    cl_uint num_platforms = 0;

    clStatus = clGetPlatformIDs(0, nullptr, &num_platforms);
    CL_ERR_FAIL_COND_MSG(clStatus, "No platforms found.");
    ERR_FAIL_COND_MSG(num_platforms == 0, "No platforms found.");

    char platform_vendor[256];
    char device_version[256];
    char lang_version[256];
    
    cl_context ctx = nullptr;
    cl_command_queue queue = nullptr;

    std::vector<cl_platform_id> platform_ids(num_platforms);
    clGetPlatformIDs(num_platforms, platform_ids.data(), nullptr);

    float best_device_score = 0;
    cl_platform_id selected_platform = nullptr;
    cl_device_id selected_device = nullptr;

    for (int i = 0; i < num_platforms; i++)
    {
        float score = 1.0f;
        cl_platform_id platform_id = platform_ids[i];

        clStatus = clGetPlatformInfo(platform_id, CL_PLATFORM_VENDOR, sizeof(platform_vendor), platform_vendor, nullptr);
        CL_ERR_FAIL_COND_MSG(clStatus, "Could not get platform info.");

        bool is_intel = strncmp("Intel", platform_vendor, strlen("Intel")) == 0;
        score *= is_intel ? 0.1f : 1;

        cl_device_id device_id = nullptr;
        cl_device_type device_type = CL_DEVICE_TYPE_GPU;
        clStatus = clGetDeviceIDs(platform_id, device_type, 1, &device_id, nullptr);
        if (clStatus != CL_SUCCESS) {
            device_type = CL_DEVICE_TYPE_CPU;
            clStatus = clGetDeviceIDs(platform_id, device_type, 1, &device_id, nullptr);
        }

        score *= device_type == CL_DEVICE_TYPE_GPU ? 1 : 0.3f;

        if (score > best_device_score) {
            selected_device = device_id;
            best_device_score = score;
            selected_platform = platform_id;
        }
    }

    ctx = clCreateContext(0, 1, &selected_device, nullptr, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(!ctx, clStatus, "Couldn't create context.");

    queue = clCreateCommandQueueWithProperties(ctx, selected_device, 0, &clStatus);
    CL_ERR_FAIL_COND_MSG(!queue, clStatus, "Couldn't create command queue.");

    *p_deviceID = selected_device;
    *p_ctx = ctx;
    *p_queue = queue;
}

void CompileOpenCLProgram(const cl_device_id p_device_id, const cl_context p_ocl_context, const char * p_soure_file, cl_program *p_ocl_program)
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


void PrintOpenCLDeviceInfo(cl_device_id device_id, cl_context contextHdl)
{
    cl_uint		uMaxComputeUnits = 0;
    cl_uint		uMaxWorkItemDim = 0;
    size_t		uMaxWorkItemSizes[3];
    cl_uint		uMaxNumSamplers = 0;
    cl_uint		uMinBaseAddrAlignSizeBits = 0;
    cl_uint		uMinBaseAddrAlignSizeBytes = 0;
    size_t		uNumBytes = 0;
    char		pDeviceVendorString[512];
    char		pDeviceNameString[512];
    char		pDriverVersionString[512];
    char		pDeviceProfileString[512];
    char		pDeviceVersionString[512];
    char		pOpenCLCVersionString[512];
    cl_int		clStatus;

    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, sizeof(pDeviceVendorString), &pDeviceVendorString, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(pDeviceNameString), &pDeviceNameString, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    clStatus = clGetDeviceInfo(device_id, CL_DRIVER_VERSION, sizeof(pDriverVersionString), &pDriverVersionString, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_PROFILE, sizeof(pDeviceProfileString), &pDeviceProfileString, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_VERSION, sizeof(pDeviceVersionString), &pDeviceVersionString, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_OPENCL_C_VERSION, sizeof(pOpenCLCVersionString), &pOpenCLCVersionString, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &uMaxComputeUnits, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &uMaxWorkItemDim, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(uMaxWorkItemSizes), &uMaxWorkItemSizes, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    size_t	uMaxWorkGroupSize;
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &uMaxWorkGroupSize, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MEM_BASE_ADDR_ALIGN, sizeof(cl_uint), &uMinBaseAddrAlignSizeBits, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, sizeof(cl_uint), &uMinBaseAddrAlignSizeBytes, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    
    cl_uint	uMaxDeviceFrequency;
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), &uMaxDeviceFrequency, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");

    /* Doesn't work
    cl_uint	uMaxImage2DWidth;
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(cl_uint), &uMaxImage2DWidth, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    */

    cl_ulong	uLocalMemSize;
    float		fLocalMemSize;
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &uLocalMemSize, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    fLocalMemSize = (float)uLocalMemSize;

    cl_long	uMaxMemAllocSize;
    float fMaxMemAllocSize;
    clStatus = clGetDeviceInfo(device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_long), &uMaxMemAllocSize, &uNumBytes);
    CL_ERR_FAIL_COND_MSG(clStatus, "clGetDeviceInfo() query failed.");
    fMaxMemAllocSize = (float)uMaxMemAllocSize;

    std::cout << "OpenCL device info:"              << std::endl;
    std::cout << "Device vendor:                  " << pDeviceVendorString << std::endl;
    std::cout << "Device name:                    " << pDeviceNameString << std::endl;
    std::cout << "Driver version:                 " << pDriverVersionString << std::endl;
    std::cout << "Device profile:                 " << pDeviceProfileString << std::endl;
    std::cout << "Device version:                 " << pDeviceVersionString << std::endl;
    std::cout << "Device OpenCL C version:        " << pOpenCLCVersionString << std::endl;
    std::cout << "Device max compute units:       " << uMaxComputeUnits << std::endl;
    std::cout << "Device max work item dim.:      " << uMaxWorkItemDim << std::endl;
    std::cout << "Device max work item sizes:     " << uMaxWorkItemSizes[0] << ", " << uMaxWorkItemSizes[1] << ", " << uMaxWorkItemSizes[2] << std::endl;
    std::cout << "Device max work group size:     " << uMaxWorkGroupSize << std::endl;
    std::cout << "Device mem. base addr. align:   " << uMinBaseAddrAlignSizeBits << " (bits)" << std::endl;
    std::cout << "Device min datatype align size: " << uMinBaseAddrAlignSizeBytes << " (bytes)" << std::endl;
    std::cout << "Device max clock frequency:     " << uMaxDeviceFrequency << std::endl;
    //std::cout << "Device image2d max width:       " << uMaxImage2DWidth << std::endl;
    std::cout << "Device local mem. size:         " << fLocalMemSize << std::endl;
    std::cout << "Device max mem alloc size:      " << fMaxMemAllocSize << std::endl;
}

#endif // UTIL_H