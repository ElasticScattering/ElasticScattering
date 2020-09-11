#pragma once

#include "src/utils/OpenCLUtils.h"
#include "src/ElasticScattering.h"
#include <random>

#include <doctest.h>
#include "TestMacros.h"

TEST_CASE("Generic gpu/cpu precision test by performing many operations on small values")
{
	int number_of_operations = 2000;
	int buffer_size = 256;

	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(&device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "test.cl", &program);

	std::vector<double> A;
	A.clear();
	A.resize(buffer_size);

	std::uniform_real_distribution<double> unif(1e-12, 9e-12);
	std::random_device r;
	std::default_random_engine re;

	for (int i = 0; i < buffer_size; i++)
		A[i] = unif(re);

	cl_int clStatus;
	cl_mem in_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

	clStatus = clEnqueueWriteBuffer(queue, in_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

	cl_kernel main_kernel = clCreateKernel(program, "many_sqrt", &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&in_buffer);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clSetKernelArg(main_kernel, 1, sizeof(int), (void*)&number_of_operations);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clSetKernelArg(main_kernel, 2, sizeof(cl_mem), (void*)&out_buffer);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	size_t global_work_size = (size_t)buffer_size;
	size_t local_work_size = 128;

	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't start test kernel execution.");

	clStatus = clFinish(queue);

	std::vector<double> gpu_results;
	gpu_results.resize(buffer_size);
	clEnqueueReadBuffer(queue, out_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, gpu_results.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

	std::vector<double> cpu_results;
	cpu_results.resize(buffer_size);

	for (int j = 0; j < buffer_size; j++) {
		for (int i = 0; i < number_of_operations; i++) {
			cpu_results[j] += A[j] * sqrt(abs(sin((double)i)));
		}
		CHECK_ALMOST(cpu_results[j], gpu_results[j], "");
	}
}

TEST_CASE("Test double2")
{
	int buffer_size = 4096;

	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(&device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "test.cl", &program);

	size_t global_work_size = buffer_size;
	size_t local_work_size = 128;

	std::vector<double> A;
	A.clear();
	A.resize(buffer_size);

	cl_int clStatus;
	cl_mem buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

	clStatus = clEnqueueWriteBuffer(queue, buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

	cl_kernel main_kernel = clCreateKernel(program, "test_double2", &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create kernel.");

	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&buffer);
	CL_FAIL_CONDITION(clStatus, "Couldn't set argument to buffer.");

	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't start test kernel execution.");

	clStatus = clFinish(queue);

	std::vector<double> gpu_results;
	gpu_results.resize(buffer_size);
	clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, sizeof(double) * buffer_size, gpu_results.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Failed to read back result.");

	double gpu_result = 0;
	double cpu_result = 0;

	for (int j = 0; j < gpu_results.size(); j++) {
		gpu_result += gpu_results[j];

		double2 x = { (double)j, (double)j };
		x = x * 1.5;
		cpu_result += x.x + x.y;
	}

	std::cout << "[Double2] CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << abs(gpu_result - cpu_result) << std::endl;
	CHECK_ALMOST(cpu_result, gpu_result, "double2 should be the same on cpu and gpu.");
}
