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
    p_init->num_iterations = 1;
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
    
    if (argc != 4) {
        std::cout << "Usage: ElasticScattering [num iterations] [lifetime | distance | stats | conductivity] [show | no-show]" << std::endl;
        std::cout << "Usage: ElasticScattering [test]" << std::endl;
        exit(0);
    }

    p_init->num_iterations = atoi(argv[1]);

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
    //ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");


    SimulationParameters sp;
    sp.region_size        = 2e-7;
    sp.dim                = 128;
    sp.particle_count     = sp.dim * sp.dim;
    sp.particle_speed     = 7e5;
    sp.particle_mass      = 5 * M0;
    sp.impurity_count     = 50000;
    sp.impurity_radius    = 2e-9;
    sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
    sp.alpha              = PI / 4.0;
    sp.phi                = 0;// -sp.alpha - 1e-10;
    sp.magnetic_field     = 0;
    sp.angular_speed      = E * sp.magnetic_field / sp.particle_mass;
    sp.tau                = 1e-13; // 3.7e-13;
    sp.integrand_steps    = 49;
    
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

    FAIL_CONDITION(pow(sp.dim, 2) != sp.particle_count, "Particles couldn't be placed in a square grid");
    FAIL_CONDITION(sp.alpha > (PI / 4.0), "Alpha should not be greater than pi/4.");
    FAIL_CONDITION(sp.alpha <= 0, "Alpha should be positive.");
    FAIL_CONDITION(sp.angular_speed < 0, "Angular speed (w) should be positive");
    FAIL_CONDITION(sp.magnetic_field < 0, "Magnetic field strength (B) should be positive");
   
    auto es = new GPUElasticScattering();
    es->Init(init, sp);
    es->Compute();

    std::cout << "+---------------------------------------------------+" << std::endl;

    bool show_demo = false;
    bool show_another = false;
    ImVec4 im_clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        //es->Compute();
        es->Draw();

        {
            static float f = 0.0f;
            static int counter = 0;

            ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.

            ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
            ImGui::Checkbox("Demo Window", &show_demo);      // Edit bools storing our window open/close state
            ImGui::Checkbox("Another Window", &show_another);

            ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
            ImGui::ColorEdit3("clear color", (float*)&im_clear_color); // Edit 3 floats representing a color

            if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
                counter++;
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
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