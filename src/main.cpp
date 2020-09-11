#define DOCTEST_CONFIG_IMPLEMENT

#include "ElasticScattering.h"
#include "tests/test_main.h"
#include "app_main.h"

#include <string>
typedef struct
{
    bool run_tests;
} InitParameters;

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->run_tests = false;

    if (argc == 2) {
        
        std::string w;
        w.assign(argv[1], strlen(argv[1]));

        p_init->run_tests = w == "test";
    }
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);

    if (init.run_tests) test_main();
    else                app_main(argc, argv);
    
    return 0;
}