#pragma once

#include <doctest.h>
#include "tests/TestMacros.h"

#include "src/sim/Simulation.h"
#include "src/sim/Grid.h"
#include "src/SimulationConfiguration.h"

TEST_CASE("CL Sigma integration should be identical to CPU sigma integration.")
{
	auto cfg = SimulationConfiguration::ParseFromeFile("tests/data/test_images.config");
	SimulationCPU sim_cpu(cfg.particles_per_row-1, cfg.quadrant_phi_steps);
	SimulationCL  sim_cl(cfg.particles_per_row-1, cfg.quadrant_phi_steps);

	auto s = cfg.settings;
	auto grid = Grid(1245124151, s.region_size, s.region_extends, s.impurity_density, s.impurity_radius, s.max_expected_impurities_in_cell);

	sim_cpu.InitSample(grid, s, false);
	sim_cl.InitSample (grid, s, false);

	SampleMetrics sm_cpu(0, 0, 0);
	SampleMetrics sm_cl (0, 0, 0);

	auto results_cpu = sim_cpu.ComputeSigmas(12, cfg.temperatures, grid, sm_cpu);
	auto results_cl  = sim_cl .ComputeSigmas(12, cfg.temperatures, grid, sm_cl);
	REQUIRE(results_cpu.size() == results_cl.size());

	for (int i = 0; i < cfg.temperatures.size(); i++)
	{
		CHECK_ALMOST(results_cpu[i].xx, results_cl[i].xx);
		CHECK_ALMOST(results_cpu[i].xy, results_cl[i].xy)
	}
}
