#define DOCTEST_CONFIG_IMPLEMENT

#include "ElasticScattering.h"
#include "tests/test_main.h"
//#include "app_main.h"
#include "sim_main.h"

#include <string>

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->mode = ProgramMode::Interactive;
    p_init->use_gpu = true;
    p_init->dont_show_info = false;

    if (argc < 2) {
        printf("|| EScl ||");
        printf("Usage:");
        printf("\t escl [sim | app | test ]");
        exit(0);
    }
        
    std::string w;
    w.assign(argv[1], strlen(argv[1]));

    if (w == "test") {
        p_init->mode = ProgramMode::Test;
    } else if (w == "sim") {
        p_init->mode = ProgramMode::Simulation;
    } else if (w == "app") {
        p_init->mode = ProgramMode::Interactive;
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
        //case Interactive: app_main(init); break;
    }
    
    return 0;
}