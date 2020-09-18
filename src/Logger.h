#pragma once

#include "src/SimulationParameters.h"
#include <vector>

class Logger {
public:
	void LogArrays(const SimulationParameters& sp, const std::vector<double>& xs, const std::vector<double>& yss) const;

	//Logger();
};