#pragma once

#include <doctest.h>
#include "TestMacros.h"
#include "TestUtils.h"

#include "src/SimulationConfiguration.h"
#include "src/sim/Simulation.h"

#include <vector>

SampleResult GetSample(Simulation& sim, const Grid& grid, const SimulationConfiguration& cfg, const bool coherent)
{
	SampleResult sample_result(cfg.temperatures.size(), cfg.magnetic_fields.size());

	sim.InitSample(grid, cfg.settings, coherent);
	Metrics metrics;

	for (int j = 0; j < cfg.magnetic_fields.size(); j++) {
		sim.ComputeLifetimes(cfg.magnetic_fields[j], grid, metrics);

		for (int i = 0; i < cfg.temperatures.size(); i++) {
			//sample_result.results[i][j] = sim.DeriveTemperature(cfg.temperatures[i]);
			sample_result.results[i][j] = sim.DeriveTemperatureWithImages(cfg.temperatures[i]).result;
		}
	}

	return sample_result;
}

void RunTestSample(const std::string file_path, const std::string imp_path, const std::vector<DataRow>& expected_results)
{
	SimulationCPU sim(31, 7);

	auto cfg = SimulationConfiguration::ParseFromeFile(file_path);

	auto impurities = GetTestImpurities(imp_path);
	REQUIRE(impurities.size() > 0);

	auto s = cfg.settings;
	auto grid = Grid(impurities, s.region_size, s.region_extends, s.impurity_radius, s.max_expected_impurities_in_cell);

	auto results_coh = GetSample(sim, grid, cfg, true);
	auto results_inc = GetSample(sim, grid, cfg, false);

	for (int i = 0; i < cfg.temperatures.size(); i++) {
		for (int j = 0; j < cfg.magnetic_fields.size(); j++) {
			DataRow row(cfg.temperatures[i], cfg.magnetic_fields[j], results_coh.results[i][j], results_inc.results[i][j]);

			auto expected_row = expected_results[j * cfg.temperatures.size() + i];
			REQUIRE_ALMOST(row.temperature, expected_row.temperature);
			REQUIRE_ALMOST(row.magnetic_field, expected_row.magnetic_field);

			printf("(%i, %i)::\n", i, j);
			printf("coh xx: %.8f <> %.8f, r/e : %.5f\n", row.coherent.xx, expected_row.coherent.xx,     row.coherent.xx   / expected_row.coherent.xx);
			printf("coh xy: %.8f <> %.8f, r/e : %.5f\n", row.coherent.xy, expected_row.coherent.xy,     row.coherent.xy   / expected_row.coherent.xy);
			printf("inc xx: %.8f <> %.8f, r/e : %.5f\n", row.incoherent.xx, expected_row.incoherent.xx, row.incoherent.xx / expected_row.incoherent.xx);
			printf("inc xy: %.8f <> %.8f, r/e : %.5f\n", row.incoherent.xy, expected_row.incoherent.xy, row.incoherent.xy / expected_row.incoherent.xy);

			CHECK_RELATIVE(row.coherent.xx, expected_row.coherent.xx);
			CHECK_RELATIVE(row.coherent.xy, expected_row.coherent.xy);
			CHECK_RELATIVE(row.incoherent.xx, expected_row.incoherent.xx);
			CHECK_RELATIVE(row.incoherent.xy, expected_row.incoherent.xy);
		}
	}
}

TEST_CASE("Integration Test")
{
	SUBCASE("File 1")
	{
		std::vector<DataRow> expected_results{
			DataRow(0.01,  0.01000, Sigma(9.856228886642308e-03, 4.539739877342534e-05), Sigma(5.780401960788526e-03, 1.989313036095190e-05)),
			DataRow(0.10,  0.01000, Sigma(9.856228886642308e-03, 4.539739877342534e-05), Sigma(5.776887324669389e-03, 1.986051613337200e-05)),
			DataRow(1.00,  0.01000, Sigma(9.856228886642308e-03, 4.539739877342534e-05), Sigma(5.741966446865339e-03, 1.953866593381223e-05)),
			DataRow(10.00,  0.01000, Sigma(9.856228886642308e-03, 4.539739877342534e-05), Sigma(5.413969604583597e-03, 1.670624623510898e-05)),
			DataRow(100.00,  0.01000, Sigma(9.856228886642308e-03, 4.539739877342534e-05), Sigma(3.428663649969118e-03, 5.654187349068526e-06)),
			DataRow(0.01, 12.00800, Sigma(9.877742028445142e-03, 1.908588517981993e-04), Sigma(5.714508777742346e-03, 1.297378574330436e-04)),
			DataRow(0.10, 12.00800, Sigma(9.877742028445142e-03, 1.908588517981993e-04), Sigma(5.711137456001301e-03, 1.295947599281806e-04)),
			DataRow(1.00, 12.00800, Sigma(9.877742028445142e-03, 1.908588517981993e-04), Sigma(5.677629948974758e-03, 1.281760275581970e-04)),
			DataRow(10.00, 12.00800, Sigma(9.877742028445142e-03, 1.908588517981993e-04), Sigma(5.361978659882804e-03, 1.151294039386555e-04)),
			DataRow(100.00, 12.00800, Sigma(9.877742028445142e-03, 1.908588517981993e-04), Sigma(3.419596367970051e-03, 4.872394578071897e-05)),
			DataRow(0.01, 24.00600, Sigma(9.919651930160808e-03, 3.745709672260805e-04), Sigma(5.560470652952115e-03, 2.346167567213858e-04)),
			DataRow(0.10, 24.00600, Sigma(9.919651930160808e-03, 3.745709672260805e-04), Sigma(5.557374126619593e-03, 2.343600966570492e-04)),
			DataRow(1.00, 24.00600, Sigma(9.919651930160808e-03, 3.745709672260805e-04), Sigma(5.526585456171533e-03, 2.318168503195231e-04)),
			DataRow(10.00, 24.00600, Sigma(9.919651930160808e-03, 3.745709672260805e-04), Sigma(5.235427448608282e-03, 2.085330857162608e-04)),
			DataRow(100.00, 24.00600, Sigma(9.919651930160808e-03, 3.745709672260805e-04), Sigma(3.395230709304849e-03, 9.032992797378409e-05)),
			DataRow(0.01, 36.00400, Sigma(9.865707931873946e-03, 5.972742035956722e-04), Sigma(5.362480270029904e-03, 3.136560495730237e-04)),
			DataRow(0.10, 36.00400, Sigma(9.865707931873946e-03, 5.972742035956722e-04), Sigma(5.359662592645578e-03, 3.133471847615739e-04)),
			DataRow(1.00, 36.00400, Sigma(9.865707931873946e-03, 5.972742035956722e-04), Sigma(5.331636203668767e-03, 3.102821394208836e-04)),
			DataRow(10.00, 36.00400, Sigma(9.865707931873946e-03, 5.972742035956722e-04), Sigma(5.065677756540391e-03, 2.818457782296645e-04)),
			DataRow(100.00, 36.00400, Sigma(9.865707931873946e-03, 5.972742035956722e-04), Sigma(3.346310075115616e-03, 1.280357749381359e-04)),
			DataRow(0.01, 48.00200, Sigma(9.883307844137637e-03, 8.192822797326902e-04), Sigma(5.161789912586584e-03, 3.779667238245487e-04)),
			DataRow(0.10, 48.00200, Sigma(9.883307844137637e-03, 8.192822797326902e-04), Sigma(5.159203348128068e-03, 3.776080778169972e-04)),
			DataRow(1.00, 48.00200, Sigma(9.883307844137637e-03, 8.192822797326902e-04), Sigma(5.133469944390398e-03, 3.740481533440381e-04)),
			DataRow(10.00, 48.00200, Sigma(9.883307844137637e-03, 8.192822797326902e-04), Sigma(4.888740807756801e-03, 3.409431901774526e-04)),
			DataRow(100.00, 48.00200, Sigma(9.883307844137637e-03, 8.192822797326902e-04), Sigma(3.281757945504918e-03, 1.590243582219660e-04)),
			DataRow(0.01, 60.00000, Sigma(9.878904094169211e-03, 1.029389427411147e-03), Sigma(4.936675399983635e-03, 4.199591991037217e-04)),
			DataRow(0.10, 60.00000, Sigma(9.878904094169211e-03, 1.029389427411147e-03), Sigma(4.934340698300760e-03, 4.195934396104412e-04)),
			DataRow(1.00, 60.00000, Sigma(9.878904094169211e-03, 1.029389427411147e-03), Sigma(4.911105577380973e-03, 4.159601029272468e-04)),
			DataRow(10.00, 60.00000, Sigma(9.878904094169211e-03, 1.029389427411147e-03), Sigma(4.689457082137099e-03, 3.819243581139883e-04)),
			DataRow(100.00, 60.00000, Sigma(9.878904094169211e-03, 1.029389427411147e-03), Sigma(3.203487163699249e-03, 1.860648537303413e-04)) };

		RunTestSample("tests/data/test_impurities.config", "tests/data/test_impurities.dat", expected_results);
	}

	SUBCASE("File 2")
	{
		std::vector<DataRow> expected_results{
		};

		//RunTestSample("tests/data/test_impurities2.config", expected_results);
	}
}


TEST_CASE("DeriveTemperature with and without logging should return same sigma result")
{
	auto cfg = SimulationConfiguration::ParseFromeFile("tests/data/test_impurities.config");

	SimulationCPU sim(cfg.particles_per_row-1, cfg.quadrant_phi_steps);

	auto s = cfg.settings;
	auto grid = Grid(21314214, s.region_size, s.region_extends, s.impurity_density, s.impurity_radius, s.max_expected_impurities_in_cell);

	Metrics metrics;

	sim.InitSample(grid, s, true);
	sim.ComputeLifetimes(10, grid, metrics);
	auto result  = sim.DeriveTemperature(2);
	auto result2 = sim.DeriveTemperatureWithImages(2);

	printf("Factor diff XX %f\n", result.xx / result2.result.xx);
	printf("Factor diff XY %f\n", result.xy / result2.result.xy);

	CHECK_ALMOST(result.xx, result2.result.xx);
	CHECK_ALMOST(result.xy, result2.result.xy);
}
