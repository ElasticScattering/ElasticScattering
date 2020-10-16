#include "ElasticScattering.h"
#include "src/escl/common.h"

double ElasticScattering::ComputeResult(const std::vector<double>& results) {
	double total = 0;
	for (int i = 0; i < results.size(); i++)
		total += results[i];

	double z = sp.region_size / (double)(sp.dim - 2);
	double result = total * z * z / 9.0;

	if (ShouldComputeSigma(sp.mode)) {
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

void ElasticScattering::CompleteSimulationParameters(ScatteringParameters& p_sp) {
	p_sp.angular_speed = E * p_sp.magnetic_field / M;

	if (p_sp.is_incoherent == 1) {
		p_sp.tau = HBAR / (KB * p_sp.temperature);
	}
	
	double area_dim = (p_sp.region_size + p_sp.region_extends * 2);
	p_sp.impurity_count = max(1, ceil(area_dim * area_dim * p_sp.impurity_density));

	particle_count = p_sp.dim * p_sp.dim;
}



bool ElasticScattering::ImpuritySettingsChanged(const ScatteringParameters& p_sp) {
	return (sp.impurity_count != p_sp.impurity_count || sp.region_extends != p_sp.region_extends || sp.region_size != p_sp.region_size || sp.impurity_seed != p_sp.impurity_seed);
};

#ifndef NO_WINDOW
uint32_t ElasticScattering::GetTextureID() const {
	return ogl.tex;
}
#endif
