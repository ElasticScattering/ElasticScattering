#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
#include "Constants.h"
#include "common_structs.h"

enum class Mode {
	AVG_LIFETIME,
	SIGMA_XX
};

typedef struct
{
	bool run_tests;
	bool show_info;
	Mode mode;
} InitParameters;

class ElasticScattering {
protected:
	std::vector<v2> impurities;
	std::vector<double> main_buffer;
	std::vector<float> pixels;

	SimulationParameters *sp = nullptr;
	SimulationParameters *last_sp;
	Mode mode;

public:
	virtual void Init(bool show_info = false) = 0;
	virtual double Compute(Mode p_mode, const SimulationParameters* p_sp) = 0;
	virtual std::vector<float> GetPixels() { return pixels; };
};

class CPUElasticScattering : public ElasticScattering {
	void PrepareCompute(const SimulationParameters* p_sp);

	double ComputeA(const v2 pos, const v2 vel);
	double ComputeB(const v2 pos, const v2 vel);

	void MakeTexture(const SimulationParameters sp)
	{
		double itau = 1.0 / sp.particle_max_lifetime; //sp.particle_count; 
		pixels.clear();
		pixels.resize(sp.particle_count * 3L);
		size_t j = 0;
		for (int i = 0; i < sp.particle_count; i++)
		{
			float k = float(main_buffer[i] * itau);
			/*

			if (k == 0) {
				pixels[j] = 1.0f;
				pixels[j + 1L] = 138.0 / 255.0;
				pixels[j + 2L] = 1.0 / 255.0;
			}
			*/
			if (k == 0) {
				pixels[j] = 0.0f;
				pixels[j + 1L] = 0.0f;
				pixels[j + 2L] = 0.0f;
			}
			else {
				pixels[j] = k;
				pixels[j + 1L] = k;
				pixels[j + 2L] = k;
			}
			j += 3;
		}
	}

public:
	virtual void Init(bool show_info = false);
	virtual double Compute(Mode p_mode, const SimulationParameters* p_sp);
};

class GPUElasticScattering : public ElasticScattering {
	void PrepareCompute(Mode p_mode, const SimulationParameters *p_sp);

	void PrepareImpurityBuffer();
	void PrepareTexKernel();

public:
	virtual void Init(bool show_info = false);
	virtual double Compute(Mode p_mode, const SimulationParameters* p_sp);
	void Draw();

	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
