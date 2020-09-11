#pragma once

#define EPSILON 0.0001
#define EPSILON_LOW_PRECISION 0.01
#define EPSILON_HIGH 1e-20
#define CHECK_ALMOST(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON_HIGH, p_msg); }
#define CHECK_APPROX(a, b)  { CHECK(abs((a)-(b)) < EPSILON); }
#define CHECK_APPROX_MSG(a, b, p_msg)  { CHECK_MESSAGE(abs((a)-(b)) < EPSILON, p_msg); }
#define CHECK_APPROX_LOW(a, b)  { CHECK(abs((a)-(b)) < EPSILON_LOW_PRECISION); }
#define CHECK_RELATIVE(a, b) { CHECK(abs(((a)-(b)) / (b)) < EPSILON); }

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


