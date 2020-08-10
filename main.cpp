#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "ElasticScattering.h"
#include "utils/Test.h"
#include "utils/ErrorMacros.h"
#include "doctest.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>


enum class Mode {
    LIFETIME,
    STATS,
    AVG_DISTANCE,
    CONDUCTIVITY
};

typedef struct
{
    bool run_tests;
    int num_iterations;
    Mode mode;
    bool show_info;
} InitParameters;

std::string ReadShaderFile2(const char* shader_file)
{
    std::ifstream file(shader_file);
    std::stringstream sstream;
    sstream << file.rdbuf();

    std::string contents = sstream.str();
    return contents;
}

void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->num_iterations = 1;
    p_init->mode = Mode::LIFETIME;
    p_init->show_info = true;
    p_init->run_tests = false;

    if (argc == 2) {
        
        std::string w;
        w.assign(argv[1], strlen(argv[1]));

        p_init->run_tests = w == "test";
    }
    return;
    
    if (argc != 4) {
        std::cout << "Usage: ElasticScattering [num iterations] [lifetime | distance | stats | conductivity] [show | no-show]" << std::endl;
        std::cout << "Usage: ElasticScattering [test]" << std::endl;
        exit(0);
    }

    p_init->num_iterations = atoi(argv[1]);

    std::unordered_map<std::string, Mode> modes
    {
        {"lifetime",     Mode::LIFETIME},
        {"distance",     Mode::AVG_DISTANCE},
        {"stats",        Mode::STATS},
        {"conductivity", Mode::CONDUCTIVITY}
    };
    std::string key;
    key.assign(argv[2], strlen(argv[2]));
    auto iterator = modes.find(key);
    FAIL_CONDITION(iterator == modes.end(), "Couldn't understand second command line argument.");
    p_init->mode = iterator->second;

    p_init->show_info = strcmp(argv[3], "show");
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);

    if (init.run_tests)
    {
        doctest::Context context;
        context.applyCommandLine(argc, argv);
        context.setOption("no-breaks", true);
        int res = context.run();

        return 0;
    }

    GLFWwindow* window;

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    if (!glfwInit())
        return -1;

    int width = 1000;
    int height = 1000;
    window = glfwCreateWindow(width, height, "Elastic Scattering", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    RECT rect;
    GetClientRect(GetDesktopWindow(), &rect);
    glfwSetWindowPos(window, rect.right / 2 - width / 2, rect.bottom / 2 - height / 2);

    glViewport(0, 0, width, height);

    GLenum error = glewInit();
    if (error != GLEW_OK) return EXIT_FAILURE;

    SimulationParameters sp;
    sp.region_size        = 1e-6;
    sp.particle_count     = 1000'000; //100'000'000;
    sp.particle_row_count = sqrt(sp.particle_count);
    sp.particle_speed     = 7e5;
    sp.particle_mass      = 5 * M0;
    sp.impurity_count     = 100;
    sp.impurity_radius    = 1.5e-8;
    sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
    sp.alpha              = PI / 4.0;
    sp.phi                = - sp.alpha - 1e-10;
    sp.magnetic_field     = 30;
    sp.angular_speed      = E * sp.magnetic_field / sp.particle_mass;
    sp.tau                = 1e-12;
    
    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "Start region size: (" << sp.region_size << ", " << sp.region_size << ")" << std::endl;
    std::cout << "Particles:         " << sp.particle_count << std::endl;
    std::cout << "Particle speed:    " << sp.particle_speed << std::endl;
    std::cout << "Particle mass:     " << sp.particle_mass << std::endl;
    std::cout << "Impurities:        " << sp.impurity_count << std::endl;
    std::cout << "Impurity radius:   " << sp.impurity_radius << std::endl;
    std::cout << "Alpha:             " << sp.alpha << std::endl;
    std::cout << "Phi:               " << sp.phi << std::endl;
    std::cout << "Magnetic field:    " << sp.magnetic_field << std::endl;
    std::cout << "Angular speed:     " << sp.angular_speed << std::endl;
    std::cout << "Tau:               " << sp.tau << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    FAIL_CONDITION(pow(sp.particle_row_count, 2) != sp.particle_count, "Particles couldn't be placed in a square grid");
    FAIL_CONDITION(sp.alpha > (PI / 4.0), "Alpha should not be greater than pi/4.");
    FAIL_CONDITION(sp.alpha <= 0, "Alpha should be positive.");
    FAIL_CONDITION(sp.angular_speed < 0, "Angular speed (w) should be positive");
    FAIL_CONDITION(sp.magnetic_field < 0, "Magnetic field strength (B) should be positive");
   
    GPUElasticScattering* es = new GPUElasticScattering();
    es->Init(sp);
    std::cout << "+---------------------------------------------------+" << std::endl;

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glEnable(GL_TEXTURE_2D);

        es->Compute();
        es->Draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}