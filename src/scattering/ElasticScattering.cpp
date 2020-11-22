#include "ElasticScattering.h"
#include "escl/lifetime.h"
#include "escl/util.h"

double ElasticScattering::SigmaFactor(const ScatteringParameters& sp) {
	double kf      = M * sp.particle_speed / HBAR;
	double outside = (E*E * kf*kf) / (2.0 * PI*PI * M * sp.region_size*sp.region_size * C1);
	double wct     = sp.angular_speed * sp.tau;
	outside *= sp.tau / (1.0 + wct*wct);

	return outside;
};

void ElasticScattering::UpdateSimulationParameters(ScatteringParameters& sp, double magnetic_field, double temperature) {
	sp.magnetic_field = magnetic_field;
	sp.temperature = temperature;

	sp.angular_speed = E * sp.magnetic_field / M;

	if (sp.is_incoherent == 1) {
		sp.tau = HBAR / (KB * sp.temperature);
		sp.default_max_lifetime = 15.0 * sp.tau;
	}
}

void ElasticScattering::CompleteSimulationParameters(ScatteringParameters& sp) {
	sp.angular_speed = E * sp.magnetic_field / M;

	if (sp.is_incoherent == 1) {
		sp.tau = HBAR / (KB * sp.temperature);
	}
	
	{
		sp.impurity_spawn_range = { -sp.region_extends, sp.region_size + sp.region_extends };
		double area_length = sp.impurity_spawn_range.y - sp.impurity_spawn_range.x;
		sp.impurity_count = max(1, ceil(area_length * area_length * sp.impurity_density));
		sp.cells_per_row = max(ceil(sqrt(sp.impurity_count / (double)sp.max_expected_impurities_in_cell)), 1);
		sp.cell_size = (sp.impurity_spawn_range.y - sp.impurity_spawn_range.x) / (double)sp.cells_per_row;
	}
	
	{
		bool diag_regions = (sp.is_diag_regions == 1);
		bool incoherent = (sp.is_incoherent == 1);

		const double incoherent_area = sp.alpha * 2.0;
		sp.integrand_angle_area = incoherent ? incoherent_area : (PI / 2.0 - incoherent_area);
		sp.integrand_step_size = sp.integrand_angle_area / (sp.integrand_steps - 1);
		
		sp.integrand_start_angle = (incoherent ? -sp.alpha : sp.alpha);
		if (diag_regions) {
			sp.integrand_start_angle += (incoherent ? (PI / 4.0) : -(PI / 4.0));
		}
	}

	sp.default_max_lifetime = 15.0 * sp.tau;
}
