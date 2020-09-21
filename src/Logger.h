#pragma once

#include "src/SimulationParameters.h"
#include "src/SimulationResult.h"
#include <vector>

class Logger {
public:
	void LogResult(const SimulationParameters& sp, const SimulationResult& sr) const;

	//Logger();
};