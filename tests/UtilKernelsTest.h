#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/scattering/ElasticScattering.h"
#include "src/utils/OpenCLUtils.h"
#include <vector>

TEST_CASE("Sum kernel")
{
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(true, &device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "src/scattering/escl/scatter.cl", &program);

	int buffer_size = 4096;
	std::vector<double> A(buffer_size, 1);

	size_t global_work_size = buffer_size / 2;
	size_t local_work_size = 128;

	cl_int clStatus;
	cl_mem in_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

	clStatus = clEnqueueWriteBuffer(queue, in_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size / local_work_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create particle_lifetimes buffer.");

	cl_kernel main_kernel = clCreateKernel(program, "sum", &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&in_buffer);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clSetKernelArg(main_kernel, 1, sizeof(cl_mem), (void*)&out_buffer);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clSetKernelArg(main_kernel, 2, sizeof(double) * local_work_size, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't start test kernel execution.");

	clStatus = clFinish(queue);

	std::vector<double> gpu_results;
	gpu_results.resize(buffer_size);
	clEnqueueReadBuffer(queue, out_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, gpu_results.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

	double gpu_result = 0;
	for (int j = 0; j < gpu_results.size(); j++) {
		gpu_result += gpu_results[j];
	}

	double cpu_result = 0;
	for (int j = 0; j < buffer_size; j++) {
		cpu_result += A[j];
	}

	std::cout << "[Sum] CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << abs(gpu_result - cpu_result) << std::endl;
	CHECK_ALMOST(cpu_result, gpu_result, "Sums on cpu and gpu should be equal.");
}

TEST_CASE("Add weights kernel")
{
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(true, &device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "src/scattering/escl/util_kernels.cl", &program);

	cl_int clStatus;
	cl_kernel main_kernel = clCreateKernel(program, "add_simpson_weights_2d", &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");


	int dim = 8;
	int buffer_size = dim * dim;

	size_t global_work_size[2] = { (size_t)dim, (size_t)dim };
	size_t local_work_size[2] = { 8, 8 };

	std::vector<double> A;
	A.clear();
	A.resize(buffer_size);

	for (int i = 0; i < buffer_size; i++)
		A[i] = 1.0;

	cl_mem buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

	clStatus = clEnqueueWriteBuffer(queue, buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&buffer);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 2, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't start test kernel execution.");

	clStatus = clFinish(queue);

	std::vector<double> gpu_results;
	gpu_results.resize(buffer_size);
	clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, sizeof(double) * buffer_size, gpu_results.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

	double gpu_result, cpu_result;
	double total_cpu = 0.0, total_gpu = 0.0;
	int limit = dim - 1;

	for (int j = 0; j < limit; j++) {
		for (int i = 0; i < limit; i++) {
			cpu_result = A[j * dim + i] * SimpsonWeight2D(i, j, limit);
			gpu_result = gpu_results[j * dim + i];

			CHECK_ALMOST(cpu_result, gpu_result, "Each weight should be the same")

			total_cpu += cpu_result;
			total_gpu += gpu_result;
		}
	}

	std::cout << "[Weights] CPU: " << total_cpu << ", GPU: " << total_gpu << ", diff: " << abs(total_gpu - total_cpu) << std::endl;
	CHECK_ALMOST(total_cpu, total_gpu, "Weights on cpu and gpu should be the same.");
}
