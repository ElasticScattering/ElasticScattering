
#include "src/scattering/escl/ScatteringParameters.h"
#include "src/SimulationResult.h"

class View
{

};

class ParameterView : View
{
public:
	bool ShowView(ScatteringParameters& sp);
};

class TexturesView : View
{
public:
	void ShowView(const std::vector<uint32_t>& textures);
};

class LogView : View
{
public:
	void ShowView(const SimulationResult& sr);
};