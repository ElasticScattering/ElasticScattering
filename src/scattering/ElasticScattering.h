#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "ImpurityIndex.h"
#include "escl/ScatteringParameters.h"
#include "src/SimulationResult.h"

#include <vector>

class ElasticScattering {
protected:
	std::vector<double> raw_lifetimes;
	ScatteringParameters sp;
	
	static double SigmaFactor(const ScatteringParameters& sp);

public:
	void UpdateSimulationParameters(ScatteringParameters& sp, double temperature);
	void CompleteSimulationParameters(ScatteringParameters& p_sp);

	virtual void ComputeLifetimes(const ScatteringParameters& sp, const ImpurityIndex& grid) = 0;
	virtual IterationResult DeriveTemperature(const double temperature) = 0;
};


class ElasticScatteringCPU : public ElasticScattering {
private:
	SigmaResult ComputeSigmas(const std::vector<double>& current_lifetimes);
	std::vector<double> IntegrateParticle(const std::vector<double>& current_lifetimes);
	Sigma IntegrateResult(const std::vector<double>& current_lifetimes);

public:
	virtual void ComputeLifetimes(const ScatteringParameters& sp, const ImpurityIndex& grid) override;
	virtual IterationResult DeriveTemperature(const double temperature) override;
};



class ElasticScatteringCL : public ElasticScattering {
	void PrepareKernels(const ScatteringParameters& sp, const size_t items_in_workgroup);

public:
	virtual void ComputeLifetimes(const ScatteringParameters& sp, const ImpurityIndex& grid) override;
	void UploadImpurities(const ImpurityIndex& grid);
	virtual IterationResult DeriveTemperature(const double temperature) override;

	ElasticScatteringCL(bool use_gpu, bool show_info, int particle_count);
	~ElasticScatteringCL();
};

#endif // ELASTIC_SCATTERING_H
