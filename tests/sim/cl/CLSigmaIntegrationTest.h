#pragma once

#include <doctest.h>
#include "tests/TestMacros.h"

#include "src/sim/Simulation.h"
#include "src/sim/Grid.h"
#include "src/SimulationConfiguration.h"

TEST_CASE("CL Sigma integration should be identical to CPU sigma integration.")
{
	auto cfg = SimulationConfiguration::ParseFromeFile("tests/data/log_with_images.config");
	SimulationCPU sim_cpu(cfg.positions_per_row, cfg.particles_per_quadrant);
	SimulationCL  sim_cl(cfg.positions_per_row, cfg.particles_per_quadrant);

	auto s = cfg.settings;
	auto grid = Grid(1245124151, s.region_size, s.region_extends, s.impurity_density, s.impurity_radius, s.target_cell_population);

	sim_cpu.InitSample(grid, s, false);
	sim_cl.InitSample (grid, s, false);

	SampleMetrics sm_cpu(0, 0, 0);
	SampleMetrics sm_cl (0, 0, 0);

	auto results_cl  = sim_cl .ComputeSigmas(12, cfg.temperatures, grid, sm_cl);
	auto results_cpu = sim_cpu.ComputeSigmas(12, cfg.temperatures, grid, sm_cpu);
	REQUIRE(results_cpu.size() == results_cl.size());

	for (int i = 0; i < cfg.temperatures.size(); i++)
	{
		CHECK_RELATIVE(results_cpu[i].xx, results_cl[i].xx);
		CHECK_RELATIVE(results_cpu[i].xy, results_cl[i].xy);

		printf("Relative difference\n");
		printf("xx: %e\n", abs((results_cpu[i].xx) - (results_cl[i].xx)) / (results_cl[i].xx));
		printf("xx: %e\n", abs((results_cpu[i].xy) - (results_cl[i].xy)) / (results_cl[i].xy));

		CHECK_ALMOST(results_cpu[i].xx, results_cl[i].xx);
		CHECK_ALMOST(results_cpu[i].xy, results_cl[i].xy)
	}
}
