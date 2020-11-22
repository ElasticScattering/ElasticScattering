#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "ImpurityGridIndex.h"
#include "escl/ScatteringParameters.h"
#include "escl/v2.h"
#include "src/SimulationResult.h"

#include <vector>

class ElasticScattering {
protected:

	static double FinishSingle(ScatteringParameters& sp, std::vector<double> &buffer);
	static double SigmaFactor(ScatteringParameters& sp);

public:
	static void CompleteSimulationParameters(ScatteringParameters& p_sp);
};

class ElasticScatteringCPU : public ElasticScattering {
public:
	static SigmaResult ComputeResult(ScatteringParameters& p_sp, ImpurityGridIndex& grid);
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
