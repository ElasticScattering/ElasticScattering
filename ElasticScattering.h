#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
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

	void MakeTexture(const SimulationParameters sp)
	{
		double itau = 1.0 / sp.particle_max_lifetime; //sp.particle_count; 
		pixels.clear();
		pixels.resize(sp.particle_count * 3L);
		size_t j = 0;
		for (int i = 0; i < sp.particle_count; i++)
		{
			float k = float(lifetimes[i] * itau);
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
	virtual void Init(SimulationParameters sp) = 0;
	virtual void Compute() = 0;
	virtual std::vector<float> GetPixels() { return pixels; };
	
};

class CPUElasticScattering : public ElasticScattering {
	double ComputeA(const v2 pos, const v2 vel, const SimulationParameters sp);
	double ComputeB(const v2 pos, const v2 vel, const SimulationParameters sp);

public:
	virtual void Init(SimulationParameters sp);
	virtual void Compute();
};

class GPUElasticScattering : public ElasticScattering {
	unsigned int tex;

	void PrepareOpenCLKernels();

public:
	virtual void Init(SimulationParameters sp);
	virtual void Compute();
	void Draw();

	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
