#pragma once

#include "SimulationParameters.h"
#include "ElasticScattering.h"

class ParametersFactory {
public:
	static ScatteringParameters& GenerateSimulation();
	static ScatteringParameters& GenerateMinimal();
	static ScatteringParameters& GenerateDefault();
	static ScatteringParameters& GenerateNoImpurities();
};