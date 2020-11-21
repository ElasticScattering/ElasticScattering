#pragma once

#include "src/scattering/escl/ScatteringParameters.h"
#include <string>

class ParametersFactory {
public:
	static ScatteringParameters& GenerateSimulation();
	static ScatteringParameters& GenerateMinimal();
	static ScatteringParameters& GenerateDefault();
	static ScatteringParameters& GenerateNoImpurities();
	static ScatteringParameters& FromFile(std::string file);
};