#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "ImpurityIndex.h"
#include "escl/ScatteringParameters.h"
#include "src/SimulationResult.h"

#include <vector>

class ElasticScattering {
protected:
	static double SigmaFactor(const ScatteringParameters& sp);

public:
	void UpdateSimulationParameters(ScatteringParameters& sp, double magnetic_field, double temperature);
	void CompleteSimulationParameters(ScatteringParameters& p_sp);

	virtual IterationResult ComputeIteration(const ScatteringParameters& sp, const ImpurityIndex& grid) = 0;
};

class ElasticScatteringCPU : public ElasticScattering {
private:
	std::vector<double> ComputeLifetimes(const ScatteringParameters& sp, const ImpurityIndex& grid);
	SigmaResult ComputeSigmas(const ScatteringParameters& sp, const std::vector<double>& particle_lifetimes);
	std::vector<double> IntegrateParticle(const ScatteringParameters& sp, const std::vector<double>& particle_lifetimes);
	Sigma IntegrateResult(const ScatteringParameters& sp, const std::vector<double>& particle_lifetimes);

public:
	IterationResult ComputeIteration(const ScatteringParameters& sp, const ImpurityIndex& grid);
};

class ElasticScatteringCL : public ElasticScattering {
	void PrepareKernels(const ScatteringParameters& sp, const size_t items_in_workgroup);

public:
	void UploadImpurities(const ImpurityIndex& grid);
	IterationResult ComputeIteration(const ScatteringParameters& sp, const ImpurityIndex& grid);

	ElasticScatteringCL(bool use_gpu, bool show_info, int particle_count);
	~ElasticScatteringCL();
};

#endif // ELASTIC_SCATTERING_H
