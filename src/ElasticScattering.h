#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
#include "src/SimulationParameters.h"

#include <GL/glew.h>

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
	double FinishSigma(double res);
	double ComputeResult(const std::vector<double>& results);

	bool AnythingChanged(const SimulationParameters& p_sp);
	virtual bool PrepareCompute(SimulationParameters& p_sp) = 0;
	void GenerateImpurities(const SimulationParameters& p_sp, bool p_random = false);

	void CompleteSimulationParameters(SimulationParameters& p_sp);

public:
	virtual bool Compute(SimulationParameters &p_sp, double &result) = 0;
	void Draw();
	uint32_t GetTextureID() const;
};

class CPUElasticScattering : public ElasticScattering {
	std::vector<double> main_buffer;
	std::vector<float>  pixels;

	virtual bool PrepareCompute(SimulationParameters &p_sp) override;
	void MakeTexture();

public:
	virtual bool Compute(SimulationParameters &p_sp, double& result) override;

	CPUElasticScattering();
};

class GPUElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(SimulationParameters &p_sp) override;
	void PrepareTexKernel(int pixels);

public:
	virtual bool Compute(SimulationParameters &p_sp, double& result) override;

	GPUElasticScattering(bool show_info = false);
	~GPUElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
