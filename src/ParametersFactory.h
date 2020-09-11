#pragma once

#include "ElasticScattering.h"

class ParametersFactory {
public:
	static SimulationParameters GenerateMinimal();
	static SimulationParameters GenerateDefault();
	static SimulationParameters GenerateNoImpurities();
};