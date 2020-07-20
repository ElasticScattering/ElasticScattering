#include "util.h"
#include <vector>
#include <fstream>
#include <sstream>

#define PROGRAM_FILE "program.cl"

cl_platform_id selected_platform;
cl_device_id selected_device;

void select_best_device();

int main(void)
{
    print_all_devices();
    cl_int clStatus;
    
    select_best_device();
    ERR_FAIL_COND(!selected_device);
    
    cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)selected_platform, 0 };
    cl_context ctx = clCreateContext(properties, 1, &selected_device, nullptr, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(!ctx, clStatus, "Couldn't create context.");

    cl_command_queue queue = clCreateCommandQueueWithProperties(ctx, selected_device, 0, &clStatus);
    CL_ERR_FAIL_COND_MSG(!queue, clStatus, "Couldn't create command queue.");

    std::fstream kernelFile(PROGRAM_FILE);
    std::string content((std::istreambuf_iterator<char>(kernelFile)), std::istreambuf_iterator<char>());
    const char* code = new char[content.size()];
    code = content.c_str();

    cl_program program = clCreateProgramWithSource(ctx, 1, (const char **) &code, nullptr, &clStatus);
    CL_ERR_FAIL_COND_MSG(!program, clStatus, "Couldn't create program.");

    clStatus = clBuildProgram(program, 1, &selected_device, NULL, NULL, NULL);
    if (clStatus != CL_SUCCESS) {
        size_t len;
        char buffer[2048];

        clGetProgramBuildInfo(program, selected_device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        CL_ERR_FAIL_COND_MSG(true, clStatus, buffer);
    }

    cl_kernel kernel1 = clCreateKernel(program, "vecadd", &clStatus);
    CL_ERR_FAIL_COND_MSG(!kernel1, clStatus, "Couldn't create kernel.");

    // Buffers, etc.


    // Free resources
    clReleaseKernel(kernel1);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(ctx);

    return 0;
}

void select_best_device() {
    cl_int clStatus;
    cl_uint num_platforms;
    size_t info_size;

    clStatus = clGetPlatformIDs(0, nullptr, &num_platforms);
    std::vector<cl_platform_id> platform_ids(num_platforms);
    clGetPlatformIDs(num_platforms, platform_ids.data(), nullptr);

    float best_device_score = 0;
    for (int i = 0; i < num_platforms; i++)
    {
        auto p = platform_ids[i];

        float score = 1.0f;

        {
            auto p_info = clGetPlatformInfo(p, CL_PLATFORM_NAME, 0, nullptr, &info_size);
            char* info = (char*)_malloca(sizeof(char) * info_size);
            clGetPlatformInfo(p, CL_PLATFORM_NAME, info_size, info, nullptr);
            bool is_intel = strncmp("Intel", info, strlen("Intel")) == 0;
            score *= is_intel ? 0.1f : 1;
        }

        cl_uint device_count = 0;
        
        cl_device_type device_type = CL_DEVICE_TYPE_GPU;
        clGetDeviceIDs(p, device_type, 0, nullptr, &device_count);
        if (device_count == 0) {
            device_type = CL_DEVICE_TYPE_CPU;
            clGetDeviceIDs(p, device_type, 0, nullptr, &device_count);
        }

        std::vector<cl_device_id> device_ids(device_count);
        clGetDeviceIDs(p, device_type, device_count, device_ids.data(), nullptr);

        score *= device_type == CL_DEVICE_TYPE_GPU ? 1 : 0.3f;

        if (score > best_device_score) {
            selected_device = device_ids[0];
            best_device_score = score;
            selected_platform = p;
        }
    }
}
