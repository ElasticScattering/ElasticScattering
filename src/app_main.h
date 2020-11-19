#pragma once

int app_main(const InitParameters &init);
void SetupContexts();
void ShutDown();
void ProcessInput();
void MainLoop(ElasticScattering& es, ScatteringParameters& sp);