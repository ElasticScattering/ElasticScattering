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

	SimulationParameters *sp = nullptr;
	SimulationParameters *last_sp;

	double FinishSigmaXX(double res) {
		double kf = sp->particle_mass * sp->particle_speed / HBAR;
		double outside = E * E * kf * kf / (2.0 * PI2 * sp->particle_mass * sp->region_size * sp->region_size * C1);
		double v = E * sp->magnetic_field * sp->tau / sp->particle_mass;
		outside /= (1.0 + v * v);

		return outside * res;
	};

	void GenerateImpurities() {
		impurities.clear();
		impurities.resize(sp->impurity_count);

		std::uniform_real_distribution<double> unif(-sp->region_extends, sp->region_size + sp->region_extends);
		std::random_device r;
		std::default_random_engine re(0);

		for (int i = 0; i < sp->impurity_count; i++)
			impurities[i] = { unif(re), unif(re) };
	};

	double ComputeResult(const std::vector<double> &results) {
		double total = 0;
		for (int i = 0; i < results.size(); i++)
			total += results[i];

		double z = sp->region_size / (sp->dim - 2);
		double result = total * z * z / 9.0;

		if (IsSigma(sp->mode))
			result = FinishSigmaXX(result);
		else
			result /= pow(sp->region_size, 2.0);

		return result * sp->particle_speed;
	}

public:
	virtual void Init(bool show_info = false) = 0;
	virtual double Compute(const SimulationParameters* p_sp) = 0;
};

class CPUElasticScattering : public ElasticScattering {
	std::vector<double> main_buffer;

	void PrepareCompute(const SimulationParameters* p_sp);

	void Lifetime();
	void LifetimePhi();

public:
	virtual void Init(bool show_info = false);
	virtual double Compute(const SimulationParameters* p_sp);
};

class GPUElasticScattering : public ElasticScattering {
	bool PrepareCompute(const SimulationParameters *p_sp);

	void PrepareImpurityBuffer();
	void PrepareTexKernel();

public:
	virtual void Init(bool show_info = false);
	virtual double Compute(const SimulationParameters* p_sp);

	void Draw();

	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
