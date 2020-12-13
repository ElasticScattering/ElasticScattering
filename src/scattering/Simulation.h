#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "Grid.h"
#include "escl/ScatteringParameters.h"
#include "src/SimulationResult.h"

#include "escl/constants.h"

#include <vector>

class Simulation {
protected:
	std::vector<double> raw_lifetimes;
	ScatteringParameters sp;
	
	double SigmaFactor(const ScatteringParameters& sp) {
		double kf = M * sp.particle_speed / HBAR;
		double outside = (E * E * kf * kf) / (2.0 * PI * PI * M * sp.region_size * sp.region_size * C1);
		double wct = sp.angular_speed * sp.tau;
		outside *= sp.tau / (1.0 + wct * wct);

		return outside;
	}

public:
	virtual void ComputeLifetimes(const ScatteringParameters& sp, const Grid& grid) = 0;
	virtual IterationResult DeriveTemperature(const double temperature) = 0;
};


class SimulationCPU : public Simulation {
private:
	SigmaResult ComputeSigmas(const std::vector<double>& current_lifetimes);
	std::vector<double> IntegrateParticle(const std::vector<double>& current_lifetimes);
	Sigma IntegrateResult(const std::vector<double>& current_lifetimes);

public:
	virtual void ComputeLifetimes(const ScatteringParameters& sp, const Grid& grid) override;
	virtual IterationResult DeriveTemperature(const double temperature) override;
};



class SimulationCL : public Simulation {
	void PrepareKernels(const ScatteringParameters& sp, const size_t items_in_workgroup);

public:
	virtual void ComputeLifetimes(const ScatteringParameters& sp, const Grid& grid) override;
	void UploadImpurities(const Grid& grid);
	virtual IterationResult DeriveTemperature(const double temperature) override;

	SimulationCL(bool use_gpu, bool show_info, int particle_count);
	~SimulationCL();
};

#endif // ELASTIC_SCATTERING_H
