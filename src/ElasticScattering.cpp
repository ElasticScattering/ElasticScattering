#include "ElasticScattering.h"

bool ElasticScattering::ImpuritySettingsChanged(const SimulationParameters &p_sp) {
	return (sp.impurity_count != p_sp.impurity_count || sp.region_extends != p_sp.region_extends || sp.region_size != p_sp.region_size || sp.impurity_seed != p_sp.impurity_seed);
};

void ElasticScattering::GenerateImpurities(const SimulationParameters &p_sp, bool p_random) {
	impurities.clear();
	impurities.resize(p_sp.impurity_count);

	std::uniform_real_distribution<double> unif(-p_sp.region_extends, p_sp.region_size + p_sp.region_extends);

	std::random_device random_device;
	unsigned int seed = p_random ? random_device() : p_sp.impurity_seed;
	std::default_random_engine re(seed);

	for (int i = 0; i < p_sp.impurity_count; i++)
		impurities[i] = { unif(re), unif(re) };
};

double ElasticScattering::FinishSigmaXX(double res) {
	double kf = M * sp.particle_speed / HBAR;
	double outside = E * E * kf * kf / (2.0 * PI2 * M * sp.region_size * sp.region_size * C1);
	double v = E * sp.magnetic_field * sp.tau / M;
	outside /= (1.0 + v * v);

	return outside * res;
};

double ElasticScattering::ComputeResult(const std::vector<double> &results) {
	double total = 0;
	for (int i = 0; i < results.size(); i++)
		total += results[i];

	double z = sp.region_size / (sp.dim - 2);
	double result = total * z * z / 9.0;

	if (IsSigma(sp.mode))
		result = FinishSigmaXX(result);
	else
		result /= pow(sp.region_size, 2.0);

	return result * sp.particle_speed;
};