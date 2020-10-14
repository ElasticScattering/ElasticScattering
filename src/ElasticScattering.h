#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
#include "datastructures/ScatteringParameters.h"
#include "datastructures/v2.h"

#ifndef NO_WINDOW
#include <GL/glew.h>
#endif

enum ProgramMode {
	Test,
	Simulation,
	Interactive
};

typedef struct
{
	ProgramMode mode;
	bool use_gpu;
	bool dont_show_info;
} InitParameters;

#ifndef NO_WINDOW
typedef struct
{
	GLuint tex;
	GLuint vbo, vao;
	GLuint shader_program;
} OpenGLResources;
#endif // NO_WINDOW

class ElasticScattering {
protected:
	std::vector<v2> impurities;

	ScatteringParameters sp;
	int particle_count;
#ifndef NO_WINDOW
	OpenGLResources ogl;
#endif
	bool first_run = true;

	bool ImpuritySettingsChanged(const ScatteringParameters& p_sp);
	double FinishSigma(double res);
	double ComputeResult(const std::vector<double>& results);

	bool AnythingChanged(const ScatteringParameters& p_sp);
	virtual bool PrepareCompute(ScatteringParameters& p_sp) = 0;
	void GenerateImpurities(const ScatteringParameters& p_sp, bool p_random = false);

	void CompleteSimulationParameters(ScatteringParameters& p_sp);

public:
	virtual bool Compute(ScatteringParameters &p_sp, double &result) = 0;
#ifndef NO_WINDOW
	uint32_t GetTextureID() const;
#endif
};

class CPUElasticScattering : public ElasticScattering {
	std::vector<double> main_buffer;
	std::vector<float>  pixels;

	virtual bool PrepareCompute(ScatteringParameters &p_sp) override;
	void MakeTexture();

public:
	virtual bool Compute(ScatteringParameters &p_sp, double& result) override;

	CPUElasticScattering();
};

class GPUElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(ScatteringParameters &p_sp) override;
	void PrepareTexKernel(int pixels);

public:
	virtual bool Compute(ScatteringParameters &p_sp, double& result) override;

	GPUElasticScattering();
	GPUElasticScattering(const InitParameters &init);
	~GPUElasticScattering();
};

class SimulationElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(ScatteringParameters& p_sp) override;

public:
	bool Compute(ScatteringParameters& p_sp, v2& result);
	virtual bool Compute(ScatteringParameters& p_sp, double& result) override;

	SimulationElasticScattering();
	SimulationElasticScattering(bool use_gpu, bool show_info);
	SimulationElasticScattering(const InitParameters& init);
	~SimulationElasticScattering();
};

#endif // ELASTIC_SCATTERING_H
