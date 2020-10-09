#define DOCTEST_CONFIG_IMPLEMENT

#include "ElasticScattering.h"
#include "tests/test_main.h"
#include "app_main.h"

#include <string>

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->mode = ProgramMode::Interactive;
    p_init->use_gpu = false;
    p_init->dont_show_info = false;

    if (argc >= 2) {
        
        std::string w;
        w.assign(argv[1], strlen(argv[1]));

        if (w == "test") {
            p_init->mode = ProgramMode::Test;
        } else if (w == "sim") {
            p_init->mode = ProgramMode::Simulation;
        } else {
            p_init->mode = ProgramMode::Interactive;
        }

        if (argc >= 3) {
            w.assign(argv[2], strlen(argv[2]));
            p_init->dont_show_info = (w == "silent");
            p_init->use_gpu = (w == "gpu"); // todo...
        }
    }
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);

    if (init.mode == ProgramMode::Test) test_main();
    if (init.mode == ProgramMode::Simulation) sim_main();
    else                app_main(init);
    
    return 0;
}