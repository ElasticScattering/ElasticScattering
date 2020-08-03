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
	std::vector<cl_double2> impurities;
	std::vector<double> lifetimes;
	std::vector<float> pixels;

	SimulationParameters sp;

public:
	virtual void Init(SimulationParameters sp) = 0;
	virtual void Compute() = 0;
	virtual std::vector<float> GetPixels() { return pixels; };
};

class CPUElasticScattering : public ElasticScattering {
	double ComputeA(const cl_double2 pos, const cl_double2 vel, const SimulationParameters sp);
	double ComputeB(const cl_double2 pos, const cl_double2 vel, const SimulationParameters sp);
	double GetBoundTime(const bool is_electron, const bool is_future) const;
	cl_double2 GetCyclotronOrbit(const cl_double2 p, const cl_double2 velocity, const double radius, const double vf, const bool is_electron) const;
	bool CirclesCross(const cl_double2 p1, const double r1, const cl_double2 p2, const double r2) const;
	std::pair<cl_double2, cl_double2> GetCrossPoints(const cl_double2 p, const double radius, const cl_double2 p2, const double r2) const;
	double GetCrossTime(const cl_double2 center, const cl_double2 pos, const cl_double2 ip, const double r, const double ir, const double w, const double clockwise) const;
	double GetPhi(const cl_double2 pos, const cl_double2 center, const double radius) const;
	double GetCrossAngle(const double p, const double q, const bool clockwise) const;

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

	void PrepareOpenCLKernels(std::vector<cl_double2> impurities, int particle_count);

public:
	virtual void Init(SimulationParameters sp);
	virtual void Compute();

	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
