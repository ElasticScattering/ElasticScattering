#ifndef UTIL_H
#define UTIL_H

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

std::string cl_error_string(int err) {
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

#define ERR_FAIL_COND_MSG(m_cond, m_msg)                                                                                 \
	if (m_cond) {                                                                                              \
		_err_print_error(FUNCTION_STR, __FILE__, __LINE__, "Condition \"" _STR(m_cond) "\" is true.", DEBUG_STR(m_msg)); \
		return;                                                                                                          \
	} else                                                                                                               \
		((void)0)

#define CL_ERR_FAIL_COND_MSG(cond, err, msg)                                                                                                         \
if (cond) {                                                                                                                                          \
    std::cout << "\x1B[31mError [" << cl_error_string(err) << "]\x1B[0m: " << msg << "\n\tat: " << "L" << __LINE__ << ": " << __FILE__ << std::endl; \
    return -1;                                                                                                                                       \
} else                                                                                                                                               \
    ((void)0);


#define ERR_FAIL_COND(cond)                                                                               \
	if (cond) {                                                                                           \
		std::cout << "\x1B[31mError\x1B[0m: \n\tat: " << "L"<< __LINE__ << ": " << __FILE__ << std::endl; \
		return -1;                                                                                        \
	} else                                                                                                \
		((void)0)

#define ERR_FAIL_COND_MSG(cond, msg)                                                                                   \
	if (cond) {                                                                                                        \
		std::cout << "\x1B[31mError\x1B[0m: " << msg << "\n\tat: " << "L"<< __LINE__ << ": " << __FILE__ << std::endl; \
		return -1;                                                                                                     \
	} else                                                                                                             \
		((void)0)

const std::vector<cl_platform_info> platform_attributes = {
    CL_PLATFORM_NAME,
    CL_PLATFORM_VENDOR,
    CL_PLATFORM_VERSION,
    CL_PLATFORM_PROFILE,
    CL_PLATFORM_EXTENSIONS
};

const std::vector<std::string> platform_attribute_names = {
    "Name",
    "Vendor",
    "Version",
    "Profile",
    "Extensions"
};

const std::vector<cl_device_info> device_attributes = {
    CL_DEVICE_NAME,
    CL_DEVICE_VENDOR,
    CL_DRIVER_VERSION,
    CL_DEVICE_VERSION,
};

const std::vector<std::string> device_attribute_names = {
    "Name",
    "Vendor",
    "Driver Version",
    "Device Version",
};

void print_all_devices() {
    cl_int clStatus;
    cl_uint num_platforms;
    size_t info_size;

    clStatus = clGetPlatformIDs(0, nullptr, &num_platforms);
    std::vector<cl_platform_id> platform_ids(num_platforms);
    clGetPlatformIDs(num_platforms, platform_ids.data(), nullptr);

    for (int i = 0; i < num_platforms; i++)
    {
        auto p = platform_ids[i];
        bool is_intel = false;

        for (int j = 0; j < platform_attributes.size(); j++) {
            auto p_info = clGetPlatformInfo(p, platform_attributes[j], 0, nullptr, &info_size);
            char* info = (char*)_malloca(sizeof(char) * info_size);
            clGetPlatformInfo(p, platform_attributes[j], info_size, info, nullptr);

            std::cout << platform_attribute_names[j] << ": " << info << std::endl;

            if (platform_attributes[j] == CL_PLATFORM_NAME) {
                is_intel = strncmp("Intel", info, strlen("Intel")) == 0;
            }
        }

        cl_uint deviceIdCount = 0;
        clGetDeviceIDs(p, CL_DEVICE_TYPE_ALL, 0, nullptr, &deviceIdCount);

        std::vector<cl_device_id> device_ids(deviceIdCount);
        clGetDeviceIDs(p, CL_DEVICE_TYPE_ALL, deviceIdCount, device_ids.data(), nullptr);

        for (int j = 0; j < deviceIdCount; j++) {
            std::cout << "Device " << j << std::endl;

            for (int k = 0; k < device_attributes.size(); k++) {
                clGetDeviceInfo(device_ids[j], device_attributes[k], 0, nullptr, &info_size);
                char* info = (char*)_malloca(sizeof(char) * info_size);
                clGetDeviceInfo(device_ids[j], device_attributes[k], info_size, info, nullptr);

                std::cout << device_attribute_names[k] << ": " << info << std::endl;
            }

            int compute_units;
            clGetDeviceInfo(device_ids[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(int), &compute_units, nullptr);
            std::cout << "CL_DEVICE_MAX_COMPUTE_UNITS" << ": " << compute_units << std::endl;
        }
    }
}

#endif // UTIL_H