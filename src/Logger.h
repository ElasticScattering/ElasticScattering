#pragma once

#include "src/datastructures/SimulationConfiguration.h"
#include "src/datastructures/SimulationResult.h"
#include <vector>

class Logger {
public:
	static void LogResult(const SimulationConfiguration& sp, const SimulationResult& sr);
};