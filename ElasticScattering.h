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
		double region_size;
		int particle_count;
		int particle_row_count;
		double particle_max_lifetime;
		double particle_speed;      // v
		double particle_mass;       // m
		int impurity_count;
		double impurity_radius;     // r
		double impurity_radius_sq;  // r^2
		double temperature;

		double alpha;
		double phi;
		double magnetic_field;      // B
		double angular_speed;       // w
	} SimulationParameters;

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

	std::vector<cl_double2> impurities;
	std::vector<double> lifetimes;
	std::vector<float> pixels;

	void ParseArgs(int argc, char** argv, InitParameters* p_init);

	void CPUElasticScattering(const SimulationParameters sp, const std::vector<cl_double2> impurities, std::vector<double> &lifetimes);
	void CPUElasticScattering2(const SimulationParameters sp, const std::vector<cl_double2> impurities, std::vector<double> &lifetimes);
	double GetBoundTime(const double phi, const double w, const double alpha, const bool is_electron, const bool is_future) const;
	cl_double2 GetCyclotronOrbit(const cl_double2 p, const cl_double2 velocity, const double radius, const double vf, const bool is_electron) const;
	bool CirclesCross(const cl_double2 p1, const double r1, const cl_double2 p2, const double r2) const;
	cl_double4 GetCrossPoints(const cl_double2 p, const double radius, const cl_double2 p2, const double r2) const;
	double GetCrossTime(const cl_double2 center, const cl_double2 pos, const cl_double2 ip, const double r, const double ir, const double w, const double clockwise) const;
	double GetPhi(const cl_double2 pos, const cl_double2 center, const double radius) const;
	double GetCrossAngle(const double p, const double q, const bool clockwise) const;

	void GPUElasticScattering(size_t size);
	void PrepareOpenCLKernels(std::vector<cl_double2> impurities, int particle_count);

	void MakeTexture(const SimulationParameters sp);

public:
	void Init(int argc, char* argv[]);
	std::vector<float> GetPixels();
	void Process();
	void Cleanup();
};

#endif // ELASTIC_SCATTERING_H
