#include "ElasticScattering.h"
#include <random>

double ElasticScattering::ComputeResult(const std::vector<double>& results) {
	double total = 0;
	for (int i = 0; i < results.size(); i++)
		total += results[i];

	double z = sp.region_size / (double)(sp.dim - 2);
	double result = total * z * z / 9.0;

	if (IsSigma(sp.mode)) {
		result = FinishSigma(result);
	}
	else {
		result /= (sp.region_size * sp.region_size);
		result *= sp.particle_speed;
		result /= PI2;
	}
	
	return result;
};

double ElasticScattering::FinishSigma(double res) {
	double kf      = M * sp.particle_speed / HBAR;
	double outside = (E*E * kf*kf) / (2.0 * PI*PI * M * sp.region_size*sp.region_size * C1);
	double wct     = sp.angular_speed * sp.tau;
	outside *= sp.tau / (1.0 + wct*wct);

	return outside * res;
};

void ElasticScattering::CompleteSimulationParameters(SimulationParameters& p_sp) {
	p_sp.angular_speed = E * p_sp.magnetic_field / M;

	if (p_sp.is_incoherent == 1) {
		p_sp.tau = HBAR / (KB * p_sp.temperature);
	}
	
	double area_dim = (p_sp.region_size + p_sp.region_extends * 2);
	p_sp.impurity_count = max(1, ceil(area_dim * area_dim * p_sp.impurity_density));

	particle_count = p_sp.dim * p_sp.dim;
}

bool ElasticScattering::AnythingChanged(const SimulationParameters& p_sp) {
	bool nothing_changed = (sp.mode == p_sp.mode && sp.impurity_seed == p_sp.impurity_seed &&
		sp.region_size == p_sp.region_size && sp.dim == p_sp.dim &&
		sp.particle_speed == p_sp.particle_speed &&
		sp.impurity_count == p_sp.impurity_count && sp.impurity_radius == p_sp.impurity_radius &&
		sp.alpha == p_sp.alpha && sp.phi == p_sp.phi && sp.temperature == p_sp.temperature &&
		sp.magnetic_field == p_sp.magnetic_field && sp.tau == p_sp.tau &&
		sp.integrand_steps == p_sp.integrand_steps && sp.is_clockwise == p_sp.is_clockwise &&
		sp.region_extends == p_sp.region_extends && sp.is_clockwise == p_sp.is_clockwise
		&& sp.is_diag_regions == p_sp.is_diag_regions && sp.is_incoherent == p_sp.is_incoherent);

	return !nothing_changed;
}

bool ElasticScattering::ImpuritySettingsChanged(const SimulationParameters& p_sp) {
	return (sp.impurity_count != p_sp.impurity_count || sp.region_extends != p_sp.region_extends || sp.region_size != p_sp.region_size || sp.impurity_seed != p_sp.impurity_seed);
};

void ElasticScattering::GenerateImpurities(const SimulationParameters& p_sp, bool p_random) {
	impurities.clear();
	impurities.resize(p_sp.impurity_count);

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

uint32_t ElasticScattering::GetTextureID() const {
	return ogl.tex;
}