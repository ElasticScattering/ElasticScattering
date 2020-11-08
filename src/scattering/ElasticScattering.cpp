#include "ElasticScattering.h"
#include "escl/lifetime.h"
#include "escl/util.h"

double ElasticScattering::FinishSingle(std::vector<double> &buffer) {
	double result = 0;

	for (int i = 0; i < buffer.size(); i++)
		result += buffer[i];

	double z = sp.region_size / (double)(sp.dim - 2);
	result *= z * z / 9.0; // phi integraal, 27?

	if (ShouldComputeSigma(sp.mode)) {
		result *= SigmaFactor();
	}
	else {
		result /= (sp.region_size * sp.region_size);
		result *= sp.particle_speed;
		result /= PI2;
	}
	
	return result;
};

double ElasticScattering::SigmaFactor() const {
	double kf      = M * sp.particle_speed / HBAR;
	double outside = (E*E * kf*kf) / (2.0 * PI*PI * M * sp.region_size*sp.region_size * C1);
	double wct     = sp.angular_speed * sp.tau;
	outside *= sp.tau / (1.0 + wct*wct);

	return outside;
};

void ElasticScattering::CompleteSimulationParameters(ScatteringParameters& p_sp) {
	p_sp.angular_speed = E * p_sp.magnetic_field / M;

	if (p_sp.is_incoherent == 1) {
		p_sp.tau = HBAR / (KB * p_sp.temperature);
	}
	
	{
		p_sp.impurity_spawn_range = { -p_sp.region_extends, p_sp.region_size + p_sp.region_extends };
		double area_length = p_sp.impurity_spawn_range.y - p_sp.impurity_spawn_range.x;
		p_sp.impurity_count = max(1, ceil(area_length * area_length * p_sp.impurity_density));
		p_sp.cells_per_row = max(sqrt(sp.impurity_count / sp.max_expected_impurities_in_cell), 1);
		p_sp.cell_size = (p_sp.impurity_spawn_range.y - p_sp.impurity_spawn_range.x) / (double)p_sp.cells_per_row;
	}
	
	{
		bool diag_regions = (p_sp.is_diag_regions == 1);
		bool incoherent = (p_sp.is_incoherent == 1);

		const double incoherent_area = p_sp.alpha * 2.0;
		p_sp.integrand_angle_area = incoherent ? incoherent_area : (PI / 2.0 - incoherent_area);
		p_sp.integrand_step_size = p_sp.integrand_angle_area / (p_sp.integrand_steps - 1);
		
		p_sp.integrand_start_angle = (incoherent ? -p_sp.alpha : p_sp.alpha);
		if (diag_regions) {
			p_sp.integrand_start_angle += (incoherent ? (PI / 4.0) : -(PI / 4.0));
		}
	}

	particle_count = p_sp.dim * p_sp.dim;
}

bool ElasticScattering::ImpuritySettingsChanged(const ScatteringParameters& p_sp) {
	return (sp.impurity_count != p_sp.impurity_count || sp.region_extends != p_sp.region_extends || sp.region_size != p_sp.region_size || sp.impurity_seed != p_sp.impurity_seed);
};
