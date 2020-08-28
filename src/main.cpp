#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "ElasticScattering.h"

#include "Test.h"
#include "utils/ErrorMacros.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>



#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    
    p_init->mode = Mode::AVG_LIFETIME;
    p_init->show_info = true;
    p_init->run_tests = false;

    if (argc == 2) {
        
        std::string w;
        w.assign(argv[1], strlen(argv[1]));

        p_init->run_tests = w == "test";
        if (p_init->run_tests) {
            p_init->show_info = false;
        }
    }
    return;
    
    if (argc != 3) {
        std::cout << "Usage: ElasticScattering [lifetime | distance | stats | conductivity] [show | no-show]" << std::endl;
        std::cout << "Usage: ElasticScattering [test]" << std::endl;
        exit(0);
    }

    
    /*
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
    */

    p_init->show_info = strcmp(argv[3], "show");
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);
    
    GLFWwindow* window;

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    if (!glfwInit())
        return -1;

    int width = 1024;
    int height = 1024;
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

    if (init.run_tests)
    {
        doctest::Context context;
        context.applyCommandLine(argc, argv);
        context.setOption("no-breaks", true);
        int res = context.run();

        return 0;
    }

    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    SimulationParameters sp;
    /*
    sp.region_size        = 2e-7;
    sp.dim                = 128;
    sp.particle_speed     = 7e5;
    sp.particle_mass      = 5 * M0;
    sp.impurity_count     = 50000;
    sp.impurity_radius    = 2e-9;
    sp.alpha              = PI / 4.0;
    sp.phi                = 0;
    sp.magnetic_field     = 0;
    sp.tau                = 1e-13;
    */

    sp.integrand_steps    = 49;
    sp.clockwise          = 1; // 1 == true, 0 == false. Can't have boolean kernel arguments :(
    sp.region_size = 1e-6;
    sp.dim = 128;
    sp.particle_speed = 7e5;
    sp.particle_mass = 5 * M0;
    sp.impurity_count = 100;
    sp.impurity_radius = 1.5e-8;
    sp.alpha = PI / 4.0;
    sp.phi = 0;// -sp.alpha - 1e-10;
    sp.magnetic_field = 0;
    sp.tau = 1e-12; // 3.7e-13;
    sp.particle_count = sp.dim * sp.dim;
    sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
    sp.angular_speed = E * sp.magnetic_field / sp.particle_mass;

    auto es = new GPUElasticScattering();
    es->Init(false);
    es->Compute(Mode::AVG_LIFETIME, &sp);

    static v2      tau_bounds = { 1e-13, 1e-10 };
    static cl_int2 integrand_steps_bounds = { 3, 99 };

    static cl_int2 count_bounds = { 1, 50000 };
    static v2      radius_bounds = { 1e-10, 1e-7 };
    static v2      region_bounds = { 1e-8,  1e-4 };

    static v2      particle_speed_bounds = { 2e2, 10e8 };
    static v2      phi_bounds = { 0, PI2 };
    static v2      magnetic_field_bounds = { 0, 80 };

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {
            ImGui::Begin("Elastic scattering");

            static int m = 0;

            ImGui::RadioButton("LT", &m, (int)Mode::AVG_LIFETIME); ImGui::SameLine();
            ImGui::RadioButton("SXX", &m, (int)Mode::SIGMA_XX);
            //ImGui::RadioButton("SXY", &init.mode, Mode::SIGMA_XY); 

            ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%e", ImGuiSliderFlags_Logarithmic);
            ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%f", 1.0f);
            ImGui::SliderInt("Integrand steps", &sp.integrand_steps, integrand_steps_bounds.x, integrand_steps_bounds.y);

            ImGui::Text("Particle count");
            ImGui::RadioButton("64x64", &sp.dim, 64); ImGui::SameLine();
            ImGui::RadioButton("128x128", &sp.dim, 128); ImGui::SameLine();
            ImGui::RadioButton("256x256", &sp.dim, 256); ImGui::SameLine();
            ImGui::RadioButton("1024x1024", &sp.dim, 1024);

            ImGui::SliderScalar("Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%e", 1.0f);
            ImGui::SliderScalar("Phi", ImGuiDataType_Double, &sp.phi, &phi_bounds.x, &phi_bounds.y, "%f", 1.0f);

            ImGui::SliderInt("Count", &sp.impurity_count, count_bounds.x, count_bounds.y);
            ImGui::SliderScalar("Region", ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%e", 1.0f);
            ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%e", 1.0f);

            if (ImGui::Button("Compute")) {
                es->Compute((Mode)m, &sp);
            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

            ImGui::End();

            es->Compute((Mode)m, &sp);
            es->Draw();
        }


        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}