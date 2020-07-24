#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "OpenCLUtils.h"

class ElasticScattering {
	enum class Mode {
		LIFETIME,
		STATS,
		AVG_DISTANCE,
		CONDUCTIVITY
	};

	typedef struct
	{
		int num_iterations;
		Mode mode;
		bool show_info;
	} InitParameters;

	typedef struct
	{
		cl_device_id deviceID;
		cl_context context;
		cl_program program;
		cl_command_queue queue;
		cl_kernel kernel;

		cl_mem db;
		cl_mem impb;
		cl_mem alive_buffer;
	} OCLResources;

	OCLResources ocl;

	size_t impurity_count, particle_count;
	cl_double2* imp_data;
	bool* alive_data;

	void ParseArgs(int argc, char** argv, InitParameters* p_init);
	void GPUElasticScattering(size_t size);
	void PrepareOpenCLKernels();

public:
	void Init(int argc, char* argv[]);
	void Process();
	void Cleanup();
};

#endif // ELASTIC_SCATTERING_H
