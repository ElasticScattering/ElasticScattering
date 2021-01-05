#pragma once

#include <doctest.h>
#include "tests/TestMacros.h"
#include "src/sim/es/cell_grid.h"
#include "src/sim/es/settings.h"
#include "src/utils/OpenCLUtils.h"
#include <vector>


TEST_CASE("CL GetNextCell")
{
	int buffer_size = 8; // num intersections.

	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(true, &device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "tests/sim/cl/kernels/test_precision.cl", &program);

	Orbit orbit;

	Intersection last_intersection;

	cl_int clStatus;
	cl_mem in_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

	clStatus = clEnqueueWriteBuffer(queue, in_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create particle_lifetimes buffer.");

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

}