#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"

#include "tests/test_main.h"
#include "utils/ErrorMacros.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include "math.h"

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
    
    GLFWwindow* window;

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    if (!glfwInit())
        return -1;

    int width = 1200;
    int height = 600;
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
        test_main();
        return 0;
    }

    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    SimulationParameters sp;
    sp.integrand_steps = 25;
    sp.region_size     = 1e-5;
    sp.dim             = 128;
    sp.particle_speed  = 1.67834e5;
    sp.impurity_count  = 100;
    sp.impurity_radius = 2e-9;
    sp.alpha           = PI / 4.0;
    sp.phi             = 1.0;
    sp.magnetic_field  = 0;
    sp.tau             = 1.5e-12;
    sp.angular_speed   = E * sp.magnetic_field / M;
    sp.region_extends  = sp.particle_speed * sp.tau * 15.0;
    
    // 1 == true, 0 == false. Can't have boolean kernel arguments :(
    sp.is_diag_regions = 0;
    sp.is_clockwise    = 0; 
    sp.is_incoherent   = 1;

    sp.mode            = MODE_SIGMA_XX;
    sp.impurity_seed   = 0;

    auto es = new GPUElasticScattering();

#if 0
    sp.impurity_count = 1;
    sp.impurity_radius = 1e-16;

    for (int i = 0; i < 50; i++) {
        sp.magnetic_field = i;
        sp.mode = MODE_SIGMA_XX;
        double result = es->Compute(sp);

        sp.mode = MODE_SIGMA_XY;

        double result2 = es->Compute(sp);
        std::cout << "" << sp.magnetic_field << " " << result << " " << result2 << std::endl;
        //std::cout << "" << i << " " << 1.0 / result << " " << result2 << " " << result / (result*result + result2*result2) << std::endl;
    }
    
    return 0;
#endif

    static v2      tau_bounds            = { 1e-13, 1e-10 };
    static cl_int2 count_bounds          = { 1, 50000 };
    static v2      radius_bounds         = { 1e-9, 1e-6 };
    static v2      region_bounds         = { 1e-8,  1e-4 };
    static v2      extends_bounds        = { 1e-7,  1e-4 };

    static v2      particle_speed_bounds = { 1e6, 1e9 };
    static v2      phi_bounds            = { 0, PI2 };
    static v2      magnetic_field_bounds = { 0, 80 };
    static bool sync_immediate           = true;
    static bool is_electron              = (sp.is_clockwise == 1);
    static bool is_diag_regions = (sp.is_diag_regions == 1);
    static bool is_incoherent = (sp.is_incoherent == 1);

    double last_result = 0;
    static int imp_seed = sp.impurity_seed;

    double ress = es->Compute(sp);

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

            static int m = sp.mode;

            ImGui::RadioButton("LT", &m, (int)MODE_DIR_LIFETIME); ImGui::SameLine();
            ImGui::RadioButton("PHI", &m, (int)MODE_PHI_LIFETIME); ImGui::SameLine();
            ImGui::RadioButton("S-XX", &m, (int)MODE_SIGMA_XX); ImGui::SameLine();
            ImGui::RadioButton("S-XY", &m, (int)MODE_SIGMA_XY);

            sp.mode = m;

            if (sp.mode != MODE_DIR_LIFETIME) {
                ImGui::Text("Phi steps");
                ImGui::RadioButton("7",  &sp.integrand_steps,  7); ImGui::SameLine();
                ImGui::RadioButton("13", &sp.integrand_steps, 13); ImGui::SameLine();
                ImGui::RadioButton("25", &sp.integrand_steps, 25); ImGui::SameLine();
                ImGui::RadioButton("49", &sp.integrand_steps, 49); ImGui::SameLine();
                ImGui::RadioButton("99", &sp.integrand_steps, 99);
            }
            else {
                ImGui::SliderScalar("Phi", ImGuiDataType_Double, &sp.phi, &phi_bounds.x, &phi_bounds.y, "%.2f", 1.0f);
            }

            ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%.2e");
            ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%.2f");
            ImGui::Checkbox("Clockwise", &is_electron);
            sp.is_clockwise = is_electron ? 1 : 0;

            if (sp.mode != MODE_DIR_LIFETIME) {
                ImGui::Checkbox("Diag. regions", &is_diag_regions); ImGui::SameLine();
                ImGui::Checkbox("Incoherent", &is_incoherent);

                sp.is_diag_regions = is_diag_regions ? 1 : 0;
                sp.is_incoherent = is_incoherent ? 1 : 0;
            }

            ImGui::Dummy(ImVec2(0.0f, 15.0f));

            ImGui::Text("Particle count");
            ImGui::RadioButton("64x64",     &sp.dim,   64); ImGui::SameLine();
            ImGui::RadioButton("128x128",   &sp.dim,  128); ImGui::SameLine();
            ImGui::RadioButton("256x256",   &sp.dim,  256); ImGui::SameLine();
            ImGui::RadioButton("1024x1024", &sp.dim, 1024);

            ImGui::SliderScalar("Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");

            ImGui::Dummy(ImVec2(0.0f, 10.0f));

            ImGui::Text("Impurities");
            ImGui::SliderInt("Count", &sp.impurity_count, count_bounds.x, count_bounds.y);
            ImGui::SliderScalar("Region", ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
            ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
            ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");
            bool impurities_updated = ImGui::Button("New seed");
            if (impurities_updated) {
                sp.impurity_seed += 1;
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            bool update = ImGui::Button("Compute"); ImGui::SameLine();
            ImGui::Checkbox("Sync immediately", &sync_immediate);

            if (sync_immediate || update || impurities_updated) {
                QueryPerformanceCounter(&beginClock);

                last_result = es->Compute(sp);
                
                //double total_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
                //std::cout << "CPU calculation time: " << total_time * 1000 << " ms" << std::endl;
                QueryPerformanceCounter(&endClock);
            }
            
            ImGui::Text("Mean free path: %.3e", last_result);
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

            ImGui::End();
        }

        es->Draw();

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