#define DOCTEST_CONFIG_IMPLEMENT

#include "SimulationConfiguration.h"
#include "SimulationRunner.h"

#ifdef _DEBUG
#include "tests/test_main.h"
#endif

#include <string>

InitParameters ParseArgs(int argc, char** argv) {
    InitParameters p_init;

    if (argc < 2) {
        p_init.config_file = "default.config";
        p_init.mode = ProgramMode::Simulation;
        return p_init;
    }

    std::string w;
    w.assign(argv[1], strlen(argv[1]));

    if (w == "help") {
        printf("Usage:\n\n");
        printf("\tElasticScattering.exe               -- Run a simulation with default.config as configuration.\n");
        printf("\tElasticScattering.exe [config_file] -- Use a custom config file.\n");
        printf("\tElasticScattering.exe help          -- Display usage.\n");
#ifdef _DEBUG
        printf("\tElasticScattering.exe test          -- Run tests.\n");
#endif
        exit(0);
    }
    if (w == "test") {
        p_init.mode = ProgramMode::Test;
    } else {
        p_init.config_file = w;
        p_init.mode = ProgramMode::Simulation;
    }

    return p_init;
}

int main(int argc, char **argv)
{
    InitParameters init = ParseArgs(argc, argv);

#ifdef _WIN32
    auto hStdin = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hStdin == INVALID_HANDLE_VALUE)
        FAIL_CONDITION(true, "GetStdHandle");

    DWORD oldConsoleMode;
    GetConsoleMode(hStdin, &oldConsoleMode);
    SetConsoleMode(hStdin, oldConsoleMode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#endif

    switch (init.mode) 
    {
        case ProgramMode::Test: 
        {
#ifdef _DEBUG
            test_main();
#endif
        } break;
        case ProgramMode::Simulation: 
        {
            SimulationRunner sr;
            sr.Run(init);
        } break;
    }
    
    return 0;
}
