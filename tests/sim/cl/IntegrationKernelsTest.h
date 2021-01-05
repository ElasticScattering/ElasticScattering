#pragma once

#include <doctest.h>
#include "tests/TestMacros.h"

#include "src/sim/es/constants.h"
#include "src/sim/es/settings.h"
#include "src/sim/es/util.h"
#include "src/utils/OpenCLUtils.h"

#include <vector>

TEST_CASE("sum kernel test")
{
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(true, &device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "src/sim/cl/integration.cl", &program);

	int buffer_size = 4096;
	double initial_value = 1.0;
	std::vector<double> A(buffer_size, initial_value);

	size_t global_work_size = buffer_size / 2;
	size_t local_work_size = 128;

	cl_int clStatus;
	cl_mem in_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size, nullptr, &clStatus);
	clEnqueueFillBuffer(queue, in_buffer, &initial_value, sizeof(double), 0, sizeof(double) * buffer_size, 0, NULL, NULL);
	//clEnqueueWriteBuffer(queue, in_buffer, CL_TRUE, 0, sizeof(double) * buffer_size, A.data(), 0, nullptr, nullptr);
	
	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size / local_work_size, nullptr, &clStatus);


	cl_kernel main_kernel = clCreateKernel(program, "sum", &clStatus);

	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&in_buffer);
	clStatus = clSetKernelArg(main_kernel, 1, sizeof(double) * local_work_size, nullptr);
	clStatus = clSetKernelArg(main_kernel, 2, sizeof(cl_mem), (void*)&out_buffer);

	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 1, nullptr, &global_work_size, &local_work_size, 0, nullptr, nullptr);
	clStatus = clFinish(queue);

	std::vector<double> gpu_results(global_work_size/local_work_size);
	clEnqueueReadBuffer(queue, out_buffer, CL_TRUE, 0, sizeof(double) * gpu_results.size(), gpu_results.data(), 0, nullptr, nullptr);

	double gpu_result = 0;
	for (int j = 0; j < gpu_results.size(); j++) {
		gpu_result += gpu_results[j];
	}

	double cpu_result = 0;
	for (int j = 0; j < buffer_size; j++) {
		cpu_result += A[j];
	}

	REQUIRE(cpu_result == buffer_size);
	CHECK_ALMOST(cpu_result, gpu_result, "Sums on cpu and gpu should be equal.");
}


TEST_CASE("3D simpson weights comparison test")
{
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(true, &device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "src/sim/cl/integration.cl", &program);

	cl_int clStatus;
	cl_kernel main_kernel = clCreateKernel(program, "apply_simpson_weights", &clStatus);
	CL_FAIL_CONDITION(clStatus, "");

	SimulationSettings ss;
	ss.positions_per_row = 7;
	ss.particles_per_quadrant = 5;
	ss.particles_per_position = ss.particles_per_quadrant * 4;
	ss.total_particles = ss.particles_per_position * ss.positions_per_row * ss.positions_per_row;

	size_t work_dim = (size_t)ss.positions_per_row + 1;
	int total_work_size = work_dim * work_dim * ss.particles_per_position;
	size_t global_work_size[3] = { work_dim, work_dim, ss.particles_per_position };
	size_t local_work_size[3] = { work_dim, 1, ss.particles_per_position };

	cl_double initial_value = 1;
	cl_mem buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * total_work_size, nullptr, &clStatus);
	clStatus = clEnqueueFillBuffer(queue, buffer, &initial_value, sizeof(double), 0, sizeof(double) * total_work_size, 0, NULL, NULL);
	
	cl_mem ss_ocl = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(SimulationSettings), nullptr, &clStatus);
	clStatus = clEnqueueWriteBuffer(queue, ss_ocl, CL_TRUE, 0, sizeof(SimulationSettings), &ss, 0, NULL, NULL);
	
	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&buffer);
	clStatus = clSetKernelArg(main_kernel, 1, sizeof(cl_mem), (void*)&ss_ocl);
	
	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
	clStatus = clFinish(queue);

	std::vector<double> gpu_results(total_work_size);
	clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, sizeof(double) * total_work_size, gpu_results.data(), 0, nullptr, nullptr);

	double cpu_total = 0;
	double gpu_total = 0;
	for (int j = 0; j < ss.positions_per_row; j++) {
		for (int i = 0; i < ss.positions_per_row; i++) {
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.particles_per_quadrant; p++) {
					cpu_total += SimpsonWeight2D(i, j, ss.positions_per_row) * SimpsonWeight(p, ss.particles_per_quadrant);
					gpu_total += gpu_results[j * work_dim * ss.particles_per_position + i * ss.particles_per_position + q * ss.particles_per_quadrant + p];
				}
			}
		}
	}
	
	CHECK_ALMOST(cpu_total, gpu_total, "CPU and CL results should be the same");
}



TEST_CASE("3D apply sigma component test")
{
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	InitializeOpenCL(true, &device, &context, &queue);

	cl_program program;
	CompileOpenCLProgram(device, context, "src/sim/cl/integration.cl", &program);

	cl_int clStatus;
	cl_kernel main_kernel = clCreateKernel(program, "apply_sigma_component", &clStatus);
	CL_FAIL_CONDITION(clStatus, "");

	SimulationSettings ss;
	ss.positions_per_row = 3;
	ss.particles_per_quadrant = 5;
	ss.particles_per_position = ss.particles_per_quadrant * 4;
	ss.total_particles = ss.particles_per_position * ss.positions_per_row * ss.positions_per_row;
	ss.signed_angular_speed = 1e7;

	ParticleSettings ps;
	ps.phi_start = 0;
	ps.phi_step_size = 0.121;

	size_t work_dim = (size_t)ss.positions_per_row + 1;
	int total_work_size = work_dim * work_dim * ss.particles_per_position;
	size_t global_work_size[3] = { work_dim, work_dim, ss.particles_per_position };
	size_t local_work_size[3] = { work_dim, 1, ss.particles_per_position };

	cl_double initial_value = 1e-14;
	cl_mem buffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(double) * total_work_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clEnqueueFillBuffer(queue, buffer, &initial_value, sizeof(double), 0, sizeof(double) * total_work_size, 0, NULL, NULL);
	CL_FAIL_CONDITION(clStatus, "");

	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(double) * total_work_size, nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "");

	cl_mem ss_ocl = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(SimulationSettings), nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clEnqueueWriteBuffer(queue, ss_ocl, CL_TRUE, 0, sizeof(SimulationSettings), &ss, 0, NULL, NULL);
	CL_FAIL_CONDITION(clStatus, "");

	cl_mem ps_ocl = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ParticleSettings), nullptr, &clStatus);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clEnqueueWriteBuffer(queue, ps_ocl, CL_TRUE, 0, sizeof(ParticleSettings), &ps, 0, NULL, NULL);
	CL_FAIL_CONDITION(clStatus, "");

	double tau = 1.0;
	int mode = MODE_SIGMA_XX;

	clStatus = clSetKernelArg(main_kernel, 0, sizeof(cl_mem), (void*)&buffer);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clSetKernelArg(main_kernel, 1, sizeof(cl_mem), (void*)&ss_ocl);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clSetKernelArg(main_kernel, 2, sizeof(cl_mem), (void*)&ps_ocl);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clSetKernelArg(main_kernel, 3, sizeof(double), (void*)&tau);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clSetKernelArg(main_kernel, 4, sizeof(int),    (void*)&mode);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clSetKernelArg(main_kernel, 5, sizeof(cl_mem), (void*)&out_buffer);
	CL_FAIL_CONDITION(clStatus, "");

	clStatus = clEnqueueNDRangeKernel(queue, main_kernel, 3, nullptr, global_work_size, local_work_size, 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "");
	clStatus = clFinish(queue);

	std::vector<double> gpu_results(total_work_size);
	clEnqueueReadBuffer(queue, out_buffer, CL_TRUE, 0, sizeof(double) * total_work_size, gpu_results.data(), 0, nullptr, nullptr);
	CL_FAIL_CONDITION(clStatus, "");

	double cpu_total = 0;
	double gpu_total = 0;
	for (int j = 0; j < ss.positions_per_row; j++) {
		for (int i = 0; i < ss.positions_per_row; i++) {
			for (int q = 0; q < 4; q++) {
				for (int p = 0; p < ss.particles_per_quadrant; p++) {
					double phi = ps.phi_start + q * HALF_PI + p * ps.phi_step_size;
					double f = (mode == MODE_SIGMA_XX) ? cos(phi) : sin(phi);
					cpu_total += f * GetSigma(initial_value, phi, tau, ss.signed_angular_speed);
					
					gpu_total += gpu_results[j * work_dim * ss.particles_per_position + i * ss.particles_per_position + q * ss.particles_per_quadrant + p];
				}
			}
		}
	}

	CHECK_RELATIVE(cpu_total, gpu_total, "CPU and CL results should be the same");
}
