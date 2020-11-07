#pragma once

#include "src/scattering/escl/ScatteringParameters.h"

class ParametersFactory {
public:
	static ScatteringParameters& GenerateSimulation();
	static ScatteringParameters& GenerateMinimal();
	static ScatteringParameters& GenerateDefault();
	static ScatteringParameters& GenerateNoImpurities();
};