#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "src/ElasticScattering.h"
#include "src/utils/OpenCLUtils.h"

TEST_CASE("Sum kernel")
{
	int buffer_size = 4096;

	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(&device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "util.h", &program);

	size_t global_work_size = buffer_size / 2;
	size_t local_work_size = 128;

	std::vector<double> A;
	A.clear();
	A.resize(buffer_size);

	for (int i = 0; i < buffer_size; i++)
		A[i] = 1.0;

	cl_int clStatus;
	cl_mem in_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create imp buffer.");

	clStatus = clEnqueueWriteBuffer(queue, in_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "Couldn't enqueue buffer.");

	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size / local_work_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "Couldn't create lifetimes buffer.");

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

TEST_CASE("Weights function") {
	int dim = 7;

	double w = GetWeight1D(0, dim);
	CHECK(w == 1);

	w = GetWeight1D(dim - 1, dim);
	CHECK(w == 1);

	w = GetWeight1D(1, dim);
	CHECK(w == 4);

	w = GetWeight1D(dim - 2, dim);
	CHECK(w == 4);

	w = GetWeight1D(dim - 3, dim);
	CHECK(w == 2);


	w = GetWeight2D(dim - 1, 0, dim);
	CHECK(w == 1);

	w = GetWeight2D(0, 0, dim);
	CHECK(w == 1);

	w = GetWeight2D(1, 0, dim);
	CHECK(w == 4);

	w = GetWeight2D(2, 0, dim);
	CHECK(w == 2);

	w = GetWeight2D(3, 0, dim);
	CHECK(w == 4);

	w = GetWeight2D(0, 1, dim);
	CHECK(w == 4);

	w = GetWeight2D(1, 1, dim);
	CHECK(w == 16);

	w = GetWeight2D(3, 3, dim);
	CHECK(w == 16);

	w = GetWeight2D(1, 2, dim);
	CHECK(w == 8);

	w = GetWeight2D(2, 1, dim);
	CHECK(w == 8);

	w = GetWeight2D(0, dim - 1, dim);
	CHECK(w == 1);
}

TEST_CASE("Add weights kernel")
{
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(&device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "common.h", &program);

	cl_int clStatus;
	cl_kernel main_kernel = clCreateKernel(program, "add_integral_weights_2d", &clStatus);
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
			cpu_result = A[j * dim + i] * GetWeight2D(i, j, limit);
			gpu_result = gpu_results[j * dim + i];

			CHECK_ALMOST(cpu_result, gpu_result, "Each weight should be the same")

			total_cpu += cpu_result;
			total_gpu += gpu_result;
		}
	}

	std::cout << "[Weights] CPU: " << total_cpu << ", GPU: " << total_gpu << ", diff: " << abs(total_gpu - total_cpu) << std::endl;
	CHECK_ALMOST(total_cpu, total_gpu, "Weights on cpu and gpu should be the same.");
}

TEST_CASE("Compare to formula (no impurities)") {
	SimulationParameters sp;
	sp.dim = 128;
	sp.particle_speed = 1.68e5;
	sp.impurity_count = 0;
	sp.impurity_radius = 2e-9;
	sp.alpha = PI / 4.0;
	sp.phi = 0;
	sp.magnetic_field = 0;
	sp.tau = 1.5e-12;
	sp.integrand_steps = 9;
	sp.is_clockwise = 0;
	sp.region_size = 1e-6;
	sp.region_extends = sp.particle_speed * sp.tau;
	sp.is_diag_regions = false;
	sp.is_incoherent = true;

	sp.mode = MODE_SIGMA_XX;
	sp.impurity_seed = 0;

	auto e = new CPUElasticScattering;

	double kf = M * sp.particle_speed / HBAR;
	double n = (kf * kf) / (PI2 * C1);
	double formula = E * E * n * sp.tau / M;

	double result = e->Compute(sp);

	std::cout << "Kf: " << kf << ", n: " << n << std::endl;

	std::cout << "CPU: " << result << ", FORM: " << formula << ", diff: " << abs(formula - result) << std::endl;
	CHECK_RELATIVE(result, formula);
}