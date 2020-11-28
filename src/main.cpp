#define DOCTEST_CONFIG_IMPLEMENT

#include "SimulationConfiguration.h"
#include "tests/test_main.h"
#include "sim_main.h"

#include <string>

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->use_gpu = true;
    p_init->dont_show_info = false; // True?

    if (argc < 2) {
        p_init->config_file = "default.config";
        p_init->mode = ProgramMode::Simulation;
        return;
    }

    std::string w;
    w.assign(argv[1], strlen(argv[1]));

    if (w == "help") {
        printf("Elastic Scattering utility - Usage:\n\n");
        printf("Running a simulation:\n");
        printf("\tescl              \t--Run a simulation with default.config as configuration");
        printf("\tescl [config_file]\t--Use custom config file.");

        printf("\nRunning tests:\n");
        printf("\tescl test");
        exit(0);
    }
    if (w == "test") {
        p_init->mode = ProgramMode::Test;
    } else {
        p_init->config_file = w;
        p_init->mode = ProgramMode::Simulation;
    }
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);

    switch (init.mode) {
        case Test       : test_main();    break;
        case Simulation : sim_main(init); break;
    }
    
    return 0;
}