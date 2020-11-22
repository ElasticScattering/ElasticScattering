#define DOCTEST_CONFIG_IMPLEMENT

#include "SimulationConfiguration.h"
#include "tests/test_main.h"
#include "sim_main.h"

#include <string>

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->use_gpu = true;
    p_init->dont_show_info = false;

    /*
    if (argc < 2) {
        printf("|| EScl ||");
        printf("Usage:");
        printf("\t escl\tRun a simulation from");
        printf("\t escl [sim | test ]");
        exit(0);
    }
    */

    std::string w;
    w.assign(argv[1], strlen(argv[1]));

    if (w == "test") {
        p_init->mode = ProgramMode::Test;
    } else {
        p_init->mode = ProgramMode::Simulation;
    }

    if (argc >= 3) {
        w.assign(argv[2], strlen(argv[2]));
        p_init->dont_show_info = (w == "silent");
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