#ifndef TEST_H
#define TEST_H

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"

#include <assert.h>
#include <random>
#include <limits>

#include "doctest.h"

#define EPSILON 0.0001
#define EPSILON_LOW_PRECISION 0.01
#define EPSILON_HIGH 1e-20
#define CHECK_ALMOST(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON_HIGH, p_msg); }
#define CHECK_APPROX(a, b)  { CHECK(abs((a)-(b)) < EPSILON); }
#define CHECK_APPROX_MSG(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON, p_msg); }
#define CHECK_APPROX_LOW(a, b)  { assert(abs((a)-(b)) < EPSILON_LOW_PRECISION); }

#define CHECK_CPU_GPU_ALMOST(p_msg)                                                            \
	cpu_result = e->Compute(sp);                                            \
	gpu_result = e2->Compute(sp);                                           \
	std::cout << "CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << abs(gpu_result-cpu_result) << std::endl;  \
	CHECK_ALMOST(cpu_result, gpu_result, p_msg);  

#define CHECK_CPU_GPU_ALMOST2(p_msg)                                                            \
	cpu_result = e->Compute(sp);                                            \
	gpu_result = e2->Compute(sp);                                           \
	diff = abs(gpu_result-cpu_result) / cpu_result;								\
	std::cout << "CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << diff << std::endl;  \
	CHECK_APPROX_MSG(diff, 0, p_msg);  

#define CHECK_CPU_GPU_APPROX(p_msg)                                                            \
	cpu_result = e->Compute(sp);                                            \
	gpu_result = e2->Compute(sp);                                           \
	std::cout << "CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << abs(gpu_result-cpu_result) << std::endl;  \
	CHECK_APPROX_MSG(cpu_result, gpu_result, p_msg);  


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

	cl_mem out_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(double) * buffer_size/local_work_size, nullptr, &clStatus);
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

	std::cout << "CPU: " << cpu_result << ", GPU: " << gpu_result << ", diff: " << abs(gpu_result - cpu_result) << std::endl;
	CHECK_ALMOST(cpu_result, gpu_result, "Sums on cpu and gpu should be equal.");
}

TEST_CASE("Weights function") {
	int dim = 8;

	double w = GetWeight1D(0, dim);
	CHECK(w == 1);

	w = GetWeight1D(dim - 1, dim);
	CHECK(w == 0);

	w = GetWeight1D(1, dim);
	CHECK(w == 4);

	w = GetWeight1D(dim-2, dim);
	CHECK(w == 1);

	w = GetWeight1D(dim - 3, dim);
	CHECK(w == 4);

	w = GetWeight2D(dim-1, 0, dim);
	CHECK(w == 0);

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
	CHECK(w == 0);
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

	for (int j = 0; j < dim; j++) {
		for (int i = 0; i < dim; i++) {
			cpu_result = A[j*dim+i] * GetWeight2D(i, j, dim);
			gpu_result = gpu_results[j * dim + i];

			CHECK_ALMOST(cpu_result, gpu_result, "Each weight should be the same")

			total_cpu += cpu_result;
			total_gpu += gpu_result;
		}
	}

	std::cout << "CPU: " << total_cpu << ", GPU: " << total_gpu << ", diff: " << abs(total_gpu - total_cpu) << std::endl;
	CHECK_ALMOST(total_cpu, total_gpu, "Weights on cpu and gpu should be the same.");
}

TEST_CASE("Compare to formula (no impurities") {
	SimulationParameters sp;
	sp.dim = 64;
	sp.particle_speed = 7e5;
	sp.impurity_count = 0;
	sp.impurity_radius = 1.5e-8;
	sp.alpha = PI / 4.0;
	sp.phi = 0;
	sp.magnetic_field = 0;
	sp.tau = 1e-12;
	sp.integrand_steps = 9;
	sp.clockwise = 0;
	sp.region_size = 1e-6;
	sp.region_extends = sp.particle_speed * sp.tau;

	sp.mode = MODE_SIGMA_XX;
	sp.impurity_seed = 0;
	
	auto e = new CPUElasticScattering;

	double kf = M * sp.particle_speed / HBAR;
    double n  = (E*E * kf*kf) / (2.0 * PI*PI * C1);
    double formula = n * sp.tau / M;

	for (int i = 0; i < 50; i++) {
		sp.magnetic_field = i;
		double result = e->Compute(sp);
		
		std::cout << "CPU: " << result << ", FORM: " << formula << ", diff: " << abs(formula - result) << std::endl;
		CHECK_ALMOST(result, formula, "Result shoudl be the same.")
	}
}

TEST_CASE("Comparing kernel results on CPU and GPU")
{
	auto e  = new CPUElasticScattering();
	auto e2 = new GPUElasticScattering();
	e2->Init(false);

	SimulationParameters sp;
	sp.region_size     = 1e-6;
	sp.dim             = 64;
	sp.particle_speed  = 7e5;
	sp.impurity_count  = 100;
	sp.impurity_radius = 1.5e-8;
	sp.alpha           = PI / 4.0;
	sp.phi             = 0;
	sp.magnetic_field  = 0;
	sp.tau             = 1e-12;
	sp.integrand_steps = 9;
	sp.clockwise       = 1;
	sp.region_extends  = sp.particle_speed * sp.tau;

	sp.mode            = MODE_DIR_LIFETIME;
	sp.impurity_seed   = 0;

	double cpu_result, gpu_result, diff;

	SUBCASE("Directional lifetime") {
		CHECK_CPU_GPU_ALMOST("Default parameters")

		sp.phi = -sp.alpha - 1e-10;
		CHECK_CPU_GPU_ALMOST("Different phi")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_ALMOST("More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_ALMOST("Larger impurities")

		sp.impurity_seed = 1;
		CHECK_CPU_GPU_ALMOST("Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_ALMOST("Magnetic field on")

		sp.clockwise = 0;
		CHECK_CPU_GPU_ALMOST("Clockwise off")
	}

	SUBCASE("Phi integrated lifetime") {
		sp.mode = MODE_PHI_LIFETIME;

		CHECK_CPU_GPU_APPROX("PHI - Default parameters")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_APPROX("PHI - More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_APPROX("PHI - Larger impurities")

		sp.impurity_seed = 2;
		CHECK_CPU_GPU_ALMOST("PHI - Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_APPROX("PHI - Magnetic field on")

		sp.clockwise = 0;
		CHECK_CPU_GPU_APPROX("PHI - Clockwise off")
	}
	

	SUBCASE("Sigma XX") {
		// SIGMA //
		sp.mode = MODE_SIGMA_XX;

		CHECK_CPU_GPU_APPROX("SXX - Default parameters")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_APPROX("SXX - More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_APPROX("SXX - Larger impurities")

		sp.impurity_seed = 3;
		CHECK_CPU_GPU_ALMOST("SXX - Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_APPROX("SXX - Magnetic field on")

		sp.clockwise = 0;
		CHECK_CPU_GPU_APPROX("SXX - Clockwise off")
	}
	
	SUBCASE("Sigma XY") {
		sp.mode = MODE_SIGMA_XY;

		CHECK_CPU_GPU_APPROX("SXY - Default parameters")

		sp.impurity_count = 200;
		CHECK_CPU_GPU_APPROX("SXY - More impurities");

		sp.impurity_radius = 1.5e-7;
		CHECK_CPU_GPU_APPROX("SXY - Larger impurities")

		sp.impurity_seed = 3;
		CHECK_CPU_GPU_ALMOST("SXY - Different impurity seed")

		sp.magnetic_field = 30;
		CHECK_CPU_GPU_APPROX("SXY - Magnetic field on")

		sp.clockwise = 0;
		CHECK_CPU_GPU_APPROX("SXY - Clockwise off")
	}
}

TEST_CASE("Cyclotron Orbit")
{
	v2 pos = { 1e-6, 1e-6 };
	v2 vel = { 7e5, 0 };
	double wc = E * 2.5 / M0;
	double vf = sqrt(vel.x * vel.x + vel.y * vel.y);
	double radius = vf / wc;

	SUBCASE("Electron is shifted up in y with positive vx")
	{
		auto c = GetCyclotronOrbit(pos, vel, radius, vf, true);
		CHECK(c.x == pos.x);
		CHECK(c.y == pos.y + radius);
	}
	
	SUBCASE("Hole is shifted down in y with positive vx")
	{
		auto c = GetCyclotronOrbit(pos, vel, radius, vf, false);
		CHECK(c.x == pos.x);
		CHECK(c.y == pos.y - radius);
	}
}

TEST_CASE("BoundTime Direction")
{
	double phi = -0.25;
	double alpha = 1;
	double w = 0.5;

	SUBCASE("Clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, true);
		CHECK(t == 0.25);
	}

	SUBCASE("Counter-clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, true, true);
		CHECK(t == 0.75);
	}

	SUBCASE("Clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, true, false);
		CHECK(t == 0.25);
	}

	SUBCASE("Counter-clockwise")
	{
		double t = GetBoundTime(phi, w, alpha, false, false);
		CHECK(t == 0.75);
	}
}

TEST_CASE("BoundTime Sectors")
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(-0.25, w, alpha, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(PI / 2 - 0.25, w, alpha, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(PI - 0.25, w, alpha, false, true);
	CHECK(t == 0.25);

	t = GetBoundTime(3 * PI / 2 - 0.25, w, alpha, false, true);
	CHECK(t == 0.25);
}

TEST_CASE("BoundTime Corners")
{
	double alpha = 1;
	double w = 0.5;

	double t = GetBoundTime(0.5, w, alpha, false, true);
	CHECK(t == 1);

	t = GetBoundTime(0.5, w, alpha, false, false);
	CHECK(t == 0);
}

TEST_CASE("Circles Cross")
{
	CHECK(!CirclesCross({ 1, 2 }, 1,    { 5, 5 }, 1));
	CHECK(!CirclesCross({ 1, 2 }, 1,    { 1, 3 }, 5));
	CHECK( CirclesCross({ 0, 0 }, 1,    { 0, 0.5 }, 0.6));
	CHECK( CirclesCross({ 0, 0 }, 1000, { 2, 0 }, 1001));
}

TEST_CASE("Circle Crosspoints")
{
	SUBCASE("Symmetric")
	{
		auto ps = GetCrossPoints({ -1, 0 }, 1.5, { 1, 0 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };
		CHECK(p1.x == 0);
		CHECK(abs(p1.y) == sqrt(1.5 * 1.5 - 1));

		CHECK(p2.x == 0);
		CHECK(abs(p2.y) == sqrt(1.5 * 1.5 - 1));
	}
	
	SUBCASE("Somewhat Symmetric")
	{
		auto ps = GetCrossPoints({ 0, 0 }, 1, { 1, 1 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };
		CHECK_APPROX(1,			  pow(p1.x, 2)     + pow(p1.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x - 1, 2) + pow(p1.y - 1, 2));

		CHECK_APPROX(1,			  pow(p2.x, 2)	   + pow(p2.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x - 1, 2) + pow(p2.y - 1, 2));
	}

	SUBCASE("Asymmetric")
	{
		auto ps = GetCrossPoints({ 100, 0 }, 100, { -1, 0 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };

		CHECK_APPROX(pow(100, 2), pow(p1.x - 100, 2) + pow(p1.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x + 1, 2)   + pow(p1.y, 2));

		CHECK_APPROX(pow(100, 2), pow(p2.x - 100, 2) + pow(p2.y, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x + 1, 2)   + pow(p2.y, 2));
	}

	SUBCASE("Asymmetric, inside")
	{
		auto ps = GetCrossPoints({ 100, -0.5 }, 100, { 1, 0.5 }, 1.5);
		v2 p1 = { ps.x, ps.y };
		v2 p2 = { ps.z, ps.w };

		CHECK_APPROX(pow(100, 2), pow(p1.x - 100, 2) + pow(p1.y + 0.5, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p1.x - 1, 2)   + pow(p1.y - 0.5, 2));

		CHECK_APPROX(pow(100, 2), pow(p2.x - 100, 2) + pow(p2.y + 0.5, 2));
		CHECK_APPROX(pow(1.5, 2), pow(p2.x - 1, 2)   + pow(p2.y - 0.5, 2));
	}
}

TEST_CASE("Phi")
{
	const v2 p1 = { 0, 0 };

	double phi = GetPhi({ 1, 0 }, p1, 1);
	CHECK(phi == 0);

	phi = GetPhi({ 0, 1 }, p1, 1);
	CHECK(phi == PI / 2);

	phi = GetPhi({ -1, 0 }, p1, 1);
	CHECK(phi == PI);

	phi = GetPhi({ 0, -1 }, p1, 1);
	CHECK(phi == 3 * PI / 2);

	phi = GetPhi({ 0.99999, -0.0045 }, p1, 1);
	CHECK_APPROX_LOW(phi, PI2);
}

TEST_CASE("Cross Angle")
{
	double a = GetCrossAngle(6.1, 0.1, true);
	CHECK(a == 6);

	a = GetCrossAngle(6.1, 0.1, false);
	CHECK_APPROX(a, PI2 - 6.0);
}

TEST_CASE("Cross Time")
{
	const v2 center = { 0, 0 }; 
	const v2 pos    = { 3, 4 };
	const double r = 5;
	double ir = 0.1;
	double w = 2;

	double t = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, true); // @todo, pos/center omdraaien geeft GetPhi assert error!
	CHECK_APPROX_LOW(t, 0.907 / 2);

	double t2 = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, false); // @todo, pos/center omdraaien geeft GetPhi assert error!
	CHECK_APPROX(t+t2, (PI2-(ir*2.0 / r))/ w);


	ir = 0.059;
	w = 100;
	
	t = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, true);
	t2 = GetFirstCrossTime(center, pos, { 5, 0 }, r, ir, w, false);
	CHECK_APPROX(t + t2, (PI2 - (ir * 2.0 / r)) / w);
}

#endif // TEST_H