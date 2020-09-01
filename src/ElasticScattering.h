#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
#include "escl/common.h"
#include <random>

typedef struct
{
	bool run_tests;
} InitParameters;

class ElasticScattering {
protected:
	std::vector<v2> impurities;

	SimulationParameters sp;
	int particle_count;
	
	bool first_run = true;

	bool ImpuritySettingsChanged(const SimulationParameters& p_sp);
	double FinishSigmaXX(double res);
	double ComputeResult(const std::vector<double>& results);
	virtual bool PrepareCompute(const SimulationParameters& p_sp) = 0;
	void GenerateImpurities(const SimulationParameters& p_sp, bool p_random = false);

public:
	virtual double Compute(const SimulationParameters &p_sp) = 0;
};

class CPUElasticScattering : public ElasticScattering {
	std::vector<double> main_buffer;

	virtual bool PrepareCompute(const SimulationParameters &p_sp);

public:
	virtual double Compute(const SimulationParameters &p_sp);
};

class GPUElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(const SimulationParameters &p_sp);
	void PrepareTexKernel(int pixels);

public:
	void Init(bool show_info = false);
	virtual double Compute(const SimulationParameters &p_sp);
	
	void Draw();

	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
