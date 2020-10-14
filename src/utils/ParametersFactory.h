#pragma once

#include "src/datastructures/ScatteringParameters.h"

class ParametersFactory {
public:
	static ScatteringParameters& GenerateSimulation();
	static ScatteringParameters& GenerateMinimal();
	static ScatteringParameters& GenerateDefault();
	static ScatteringParameters& GenerateNoImpurities();
};