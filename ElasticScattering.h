#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "OpenCLUtils.h"
#include "Constants.h"

typedef struct
{
	double region_size;
	int particle_count;
	int particle_row_count;
	double particle_max_lifetime;
	double particle_speed;      // v
	double particle_mass;       // m
	int impurity_count;
	double impurity_radius;     // r
	double impurity_radius_sq;  // r^2
	double tau;

	double alpha;
	double phi;
	double magnetic_field;      // B
	double angular_speed;       // w
} SimulationParameters;

class ElasticScattering {
protected:
	std::vector<v2> impurities;
	std::vector<double> lifetimes;
	std::vector<float> pixels;

	SimulationParameters sp;

public:
	virtual void Init(SimulationParameters sp) = 0;
	virtual void Compute() = 0;
	virtual std::vector<float> GetPixels() { return pixels; };
};

class CPUElasticScattering : public ElasticScattering {
	double ComputeA(const v2 pos, const v2 vel, const SimulationParameters sp);
	double ComputeB(const v2 pos, const v2 vel, const SimulationParameters sp);

	void MakeTexture(const SimulationParameters sp);

public:
	virtual void Init(SimulationParameters sp);
	virtual void Compute();
};

class GPUElasticScattering : public ElasticScattering {
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

	void PrepareOpenCLKernels(int particle_count);

public:
	virtual void Init(SimulationParameters sp);
	virtual void Compute();

	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
