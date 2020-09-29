#define DOCTEST_CONFIG_IMPLEMENT

#include "ElasticScattering.h"
#include "tests/test_main.h"
#include "app_main.h"

#include <string>

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->run_tests = false;
    p_init->use_gpu = false;
    p_init->dont_show_info = false;

    if (argc >= 2) {
        
        std::string w;
        w.assign(argv[1], strlen(argv[1]));

        p_init->run_tests = (w == "test");
        p_init->use_gpu   = (w == "gpu");

        if (argc >= 3) {
            w.assign(argv[2], strlen(argv[2]));
            p_init->dont_show_info = (w == "silent");
        }
    }
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);

    if (init.run_tests) test_main();
    else                app_main(init);
    
    return 0;
}