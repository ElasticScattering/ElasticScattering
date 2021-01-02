#pragma once

#include <doctest.h>
#include "tests/TestMacros.h"
#include "tests/TestUtils.h"

#include "src/SimulationConfiguration.h"
#include "src/sim/Simulation.h"
#include "src/Logger.h"

#include <vector>

SampleResult GetSample(Simulation& sim, const Grid& grid, const SimulationConfiguration& cfg, const bool coherent, SampleMetrics& sample_metrics)
{
	SampleResult sample_result(cfg.temperatures.size(), cfg.magnetic_fields.size());
	
	sim.InitSample(grid, cfg.settings, coherent);
	for (int i = 0; i < cfg.magnetic_fields.size(); i++)
		sample_result.results[i] = sim.ComputeSigmas(cfg.magnetic_fields[i], cfg.temperatures, grid, sample_metrics);

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

	SampleMetrics metrics_coh(0, 0, 0);
	auto results_coh = GetSample(sim, grid, cfg, true, metrics_coh);
	
	SampleMetrics metrics_inc(0, 0, 0);
	auto results_inc = GetSample(sim, grid, cfg, false, metrics_inc);

	for (int i = 0; i < cfg.temperatures.size(); i++) {
		for (int j = 1; j < cfg.magnetic_fields.size(); j++) {
			DataRow row(cfg.temperatures[i], cfg.magnetic_fields[j], results_coh.results[j][i], results_inc.results[j][i]);

			auto expected_row = expected_results[j * cfg.temperatures.size() + i];
			REQUIRE_ALMOST(row.temperature, expected_row.temperature);
			REQUIRE_ALMOST(row.magnetic_field, expected_row.magnetic_field);

			/*
			printf("(%i, %i)::\n", i, j);
			printf("coh xx: %.8f <> %.8f, r/e : %.5f\n", row.coherent.xx, expected_row.coherent.xx,     row.coherent.xx   / expected_row.coherent.xx);
			printf("coh xy: %.8f <> %.8f, r/e : %.5f\n", row.coherent.xy, expected_row.coherent.xy,     row.coherent.xy   / expected_row.coherent.xy);
			printf("inc xx: %.8f <> %.8f, r/e : %.5f\n", row.incoherent.xx, expected_row.incoherent.xx, row.incoherent.xx / expected_row.incoherent.xx);
			printf("inc xy: %.8f <> %.8f, r/e : %.5f\n", row.incoherent.xy, expected_row.incoherent.xy, row.incoherent.xy / expected_row.incoherent.xy);
			*/

			CHECK_RELATIVE(row.coherent.xx, expected_row.coherent.xx);
			CHECK_RELATIVE(row.coherent.xy, expected_row.coherent.xy);
			CHECK_RELATIVE(row.incoherent.xx, expected_row.incoherent.xx);
			CHECK_RELATIVE(row.incoherent.xy, expected_row.incoherent.xy);

		}
	}

	for (int j = 0; j < cfg.magnetic_fields.size(); j++) {
		CHECK(metrics_inc.iteration_metrics[j].particles_escaped == 0);
		CHECK(metrics_coh.iteration_metrics[j].particles_escaped == 0);
	}
}

TEST_CASE("Integration Test")
{
	SUBCASE("File 1")
	{
		std::vector<DataRow> expected_results{
			DataRow(0.01,  0.01000, Sigma(9.800103732311494e-03, 7.771488577983383e-06), Sigma(5.725051164160944e-03, 2.061334235275216e-05)),
			DataRow(0.10,  0.01000, Sigma(9.800103732311494e-03, 7.771488577983383e-06), Sigma(5.721603950055542e-03, 2.057533776712238e-05)),
			DataRow(1.00,  0.01000, Sigma(9.800103732311494e-03, 7.771488577983383e-06), Sigma(5.687349018423400e-03, 2.020037489270508e-05)),
			DataRow(10.00,  0.01000, Sigma(9.800103732311494e-03, 7.771488577983383e-06), Sigma(5.365262758915634e-03, 1.690977662128921e-05)),
			DataRow(100.00,  0.01000, Sigma(9.800103732311494e-03, 7.771488577983383e-06), Sigma(3.404663308768081e-03, 4.779641758969151e-06)),
			DataRow(0.01, 12.00800, Sigma(9.826143713940532e-03, 1.811365024574828e-04), Sigma(5.678718243508259e-03, 1.357357518514256e-04)),
			DataRow(0.10, 12.00800, Sigma(9.826143713940532e-03, 1.811365024574828e-04), Sigma(5.675357976345992e-03, 1.355793079234559e-04)),
			DataRow(1.00, 12.00800, Sigma(9.826143713940532e-03, 1.811365024574828e-04), Sigma(5.641962015542396e-03, 1.340288440521277e-04)),
			DataRow(10.00, 12.00800, Sigma(9.826143713940532e-03, 1.811365024574828e-04), Sigma(5.327508695131951e-03, 1.198228259094654e-04)),
			DataRow(100.00, 12.00800, Sigma(9.826143713940532e-03, 1.811365024574828e-04), Sigma(3.397167550109725e-03, 4.932576312743141e-05)),
			DataRow(0.01, 24.00600, Sigma(9.850655229401642e-03, 3.875341528783627e-04), Sigma(5.511780042513249e-03, 2.336577813124129e-04)),
			DataRow(0.10, 24.00600, Sigma(9.850655229401642e-03, 3.875341528783627e-04), Sigma(5.508742506691057e-03, 2.334011659366217e-04)),
			DataRow(1.00, 24.00600, Sigma(9.850655229401642e-03, 3.875341528783627e-04), Sigma(5.478537145341734e-03, 2.308585434691764e-04)),
			DataRow(10.00, 24.00600, Sigma(9.850655229401642e-03, 3.875341528783627e-04), Sigma(5.192609783982774e-03, 2.075960824686504e-04)),
			DataRow(100.00, 24.00600, Sigma(9.850655229401642e-03, 3.875341528783627e-04), Sigma(3.375110665051661e-03, 8.993659417056850e-05)),
			DataRow(0.01, 36.00400, Sigma(9.824991276216539e-03, 5.939282278737960e-04), Sigma(5.313640619667328e-03, 3.089641556649546e-04)),
			DataRow(0.10, 36.00400, Sigma(9.824991276216539e-03, 5.939282278737960e-04), Sigma(5.310863450905830e-03, 3.086592170373555e-04)),
			DataRow(1.00, 36.00400, Sigma(9.824991276216539e-03, 5.939282278737960e-04), Sigma(5.283239030699984e-03, 3.056332132326723e-04)),
			DataRow(10.00, 36.00400, Sigma(9.824991276216539e-03, 5.939282278737960e-04), Sigma(5.021007850155433e-03, 2.775666566876359e-04)),
			DataRow(100.00, 36.00400, Sigma(9.824991276216539e-03, 5.939282278737960e-04), Sigma(3.321782187098294e-03, 1.261416813620387e-04)),
			DataRow(0.01, 48.00200, Sigma(9.822459176830620e-03, 8.166429449514070e-04), Sigma(5.108967240260496e-03, 3.732159133342496e-04)),
			DataRow(0.10, 48.00200, Sigma(9.822459176830620e-03, 8.166429449514070e-04), Sigma(5.106423046864812e-03, 3.728636693591780e-04)),
			DataRow(1.00, 48.00200, Sigma(9.822459176830620e-03, 8.166429449514070e-04), Sigma(5.081110525531792e-03, 3.693671742060704e-04)),
			DataRow(10.00, 48.00200, Sigma(9.822459176830620e-03, 8.166429449514070e-04), Sigma(4.840321003433200e-03, 3.368416534894883e-04)),
			DataRow(100.00, 48.00200, Sigma(9.822459176830620e-03, 8.166429449514070e-04), Sigma(3.255545305930707e-03, 1.577213278430248e-04)),
			DataRow(0.01, 60.00000, Sigma(9.799939501834578e-03, 1.021948030516029e-03), Sigma(4.895327489981128e-03, 4.189645385540096e-04)),
			DataRow(0.10, 60.00000, Sigma(9.799939501834578e-03, 1.021948030516029e-03), Sigma(4.893014418885358e-03, 4.185968655701820e-04)),
			DataRow(1.00, 60.00000, Sigma(9.799939501834578e-03, 1.021948030516029e-03), Sigma(4.869994974415824e-03, 4.149448658540921e-04)),
			DataRow(10.00, 60.00000, Sigma(9.799939501834578e-03, 1.021948030516029e-03), Sigma(4.650438063348808e-03, 3.807647647128992e-04)),
			DataRow(100.00, 60.00000, Sigma(9.799939501834578e-03, 1.021948030516029e-03), Sigma(3.179082341466999e-03, 1.851160576163891e-04))
		};

		RunTestSample("tests/data/test_smaller_case.config", "tests/data/test_impurities.dat", expected_results);
	}
}


TEST_CASE("DeriveTemperature with and without logging should return same sigma result")
{
	auto cfg = SimulationConfiguration::ParseFromeFile("tests/data/test_images.config");
	REQUIRE(cfg.temperatures.size() == 1);

	SimulationCPU sim(cfg.particles_per_row-1, cfg.quadrant_phi_steps);

	auto s = cfg.settings;
	auto grid = Grid(21314214, s.region_size, s.region_extends, s.impurity_density, s.impurity_radius, s.max_expected_impurities_in_cell);

	SampleMetrics metrics(0, 0, 0);

	sim.InitSample(grid, s, true);
	auto result  = sim.ComputeSigmas          (10, cfg.temperatures, grid, metrics)[0];
	auto result2 = sim.ComputeSigmasWithImages(10, cfg.temperatures, grid, metrics)[0].result;

	CHECK_ALMOST(result.xx, result2.xx);
	CHECK_ALMOST(result.xy, result2.xy);
}
