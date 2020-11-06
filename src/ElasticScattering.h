#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include <vector>
#include "datastructures/ScatteringParameters.h"
#include "datastructures/v2.h"
#include "datastructures/IntegralResult.h"
#include "datastructures/SimulationResult.h"
#include "ImpurityGrid.h"

#include <GL/glew.h>
#include "utils/OpenGLUtils.h"

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

typedef struct
{
	GLuint tex;
	GLuint vbo, vao;
	GLuint shader_program;
} OpenGLResources;

class ElasticScattering {
protected:
	ImpurityGrid grid;

	ScatteringParameters sp;
	int particle_count;
	
	bool first_run = true;

	virtual bool PrepareCompute(ScatteringParameters& p_sp) = 0;
	void CompleteSimulationParameters(ScatteringParameters& p_sp);
	bool ImpuritySettingsChanged(const ScatteringParameters& p_sp);

	double FinishSingle(std::vector<double> &buffer);
	double SigmaFactor() const;

public:
	virtual SigmaResult ComputeResult(ScatteringParameters& p_sp) = 0;
	virtual bool ComputeSingle(ScatteringParameters& p_sp, double &buffer) = 0;
};

class CPUElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(ScatteringParameters &p_sp) override;

public:
	virtual SigmaResult ComputeResult(ScatteringParameters &p_sp) override;
	virtual bool ComputeSingle(ScatteringParameters& p_sp, double& buffer) override;
};

/*
class GPUElasticScattering : public ElasticScattering {
	OpenGLResources ogl;

	virtual bool PrepareCompute(ScatteringParameters &p_sp) override;
	void PrepareTexKernel(int pixels);

public:
	virtual bool Compute(ScatteringParameters &p_sp, IntegralResult& result) override;

	GPUElasticScattering();
	GPUElasticScattering(const InitParameters &init);
	~GPUElasticScattering();
	
	uint32_t GetTextureID() const;
};

class SimulationElasticScattering : public ElasticScattering {
	virtual bool PrepareCompute(ScatteringParameters& p_sp) override;

public:
	bool Compute(ScatteringParameters& p_sp, v2& result);
	virtual bool Compute(ScatteringParameters& p_sp, IntegralResult& result) override;

	SimulationElasticScattering();
	SimulationElasticScattering(bool use_gpu, bool show_info);
	SimulationElasticScattering(const InitParameters& init);
	~SimulationElasticScattering();
};
*/
#endif // ELASTIC_SCATTERING_H
