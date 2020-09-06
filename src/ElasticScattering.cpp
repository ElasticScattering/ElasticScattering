#include "ElasticScattering.h"
#include <random>

double ElasticScattering::ComputeResult(const std::vector<double>& results) {
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

double ElasticScattering::FinishSigmaXX(double res) {
	double kf = M * sp.particle_speed / HBAR;
	double outside = (E*E * kf*kf) / (2.0 * PI*PI * M * sp.region_size*sp.region_size * C1);
	double v = E * sp.magnetic_field * sp.tau / M;
	outside *= sp.tau / (1.0 + v*v);

	return outside * res;
};

void ElasticScattering::CompleteSimulationParameters() {
	sp.angular_speed = E * sp.magnetic_field / M;
	particle_count = sp.dim * sp.dim;
}

bool ElasticScattering::AnythingChanged(const SimulationParameters& p_sp) {
	bool nothing_changed = (sp.mode == p_sp.mode && sp.impurity_seed == p_sp.impurity_seed &&
		sp.region_size == p_sp.region_size && sp.dim == p_sp.dim &&
		sp.particle_speed == p_sp.particle_speed &&
		sp.impurity_count == p_sp.impurity_count && sp.impurity_radius == p_sp.impurity_radius &&
		sp.alpha == p_sp.alpha && sp.phi == p_sp.phi &&
		sp.magnetic_field == p_sp.magnetic_field && sp.tau == p_sp.tau &&
		sp.integrand_steps == p_sp.integrand_steps && sp.clockwise == p_sp.clockwise &&
		sp.region_extends == p_sp.region_extends);

	return !nothing_changed;
}

bool ElasticScattering::ImpuritySettingsChanged(const SimulationParameters& p_sp) {
	return (sp.impurity_count != p_sp.impurity_count || sp.region_extends != p_sp.region_extends || sp.region_size != p_sp.region_size || sp.impurity_seed != p_sp.impurity_seed);
};

void ElasticScattering::GenerateImpurities(const SimulationParameters& p_sp, bool p_random) {
	impurities.clear();
	impurities.resize(100);

	std::uniform_real_distribution<double> unif(-p_sp.region_extends, p_sp.region_size + p_sp.region_extends);

	std::random_device random_device;
	unsigned int seed = p_random ? random_device() : p_sp.impurity_seed;
	std::default_random_engine re(seed);

	for (int i = 0; i < p_sp.impurity_count; i++)
		impurities[i] = { unif(re), unif(re) };
};

void ElasticScattering::Draw()
{
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, ogl.tex);

	glUseProgram(ogl.shader_program);
	glBindVertexArray(ogl.vao);
	glDrawArrays(GL_QUADS, 0, 4);
}