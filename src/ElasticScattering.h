#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
#include "escl/common.h"

#include <GL/glew.h>

typedef struct
{
	bool run_tests;
} InitParameters;

typedef struct
{
	GLuint tex;
	GLuint vbo, vao;
	GLuint shader_program;
} OpenGLResources;

class ElasticScattering {
protected:
	std::vector<v2> impurities;

	SimulationParameters sp;
	int particle_count;
	
	OpenGLResources ogl;
	
	bool first_run = true;

	bool ImpuritySettingsChanged(const SimulationParameters& p_sp);
	double FinishSigmaXX(double res);
	double ComputeResult(const std::vector<double>& results);
	virtual bool PrepareCompute(const SimulationParameters& p_sp) = 0;
	void GenerateImpurities(const SimulationParameters& p_sp, bool p_random = false);

	void CompleteSimulationParameters();

public:
	virtual double Compute(const SimulationParameters &p_sp) = 0;
	void Draw();
};

class CPUElasticScattering : public ElasticScattering {
	std::vector<double> main_buffer;
	std::vector<float>  pixels;

	virtual bool PrepareCompute(const SimulationParameters &p_sp);
	void MakeTexture();

public:
	virtual double Compute(const SimulationParameters &p_sp);
	
	CPUElasticScattering();
};

class GPUElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(const SimulationParameters &p_sp);
	void PrepareTexKernel(int pixels);

public:
	virtual double Compute(const SimulationParameters &p_sp);
	
	GPUElasticScattering(bool show_info = false);
	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
