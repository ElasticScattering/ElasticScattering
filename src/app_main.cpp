#pragma once
#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "app_main.h"
#include "ElasticScattering.h"
#include "utils/OpenCLUtils.h"
#include "utils/ErrorMacros.h"
#include "src/ParametersFactory.h"

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

int app_main(int argc, char** argv)
{
    GLFWwindow* window;

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    if (!glfwInit())
        return -1;

    int width = 1600;
    int height = 800;
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

    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    LARGE_INTEGER beginClock, endClock, clockFrequency;
    QueryPerformanceFrequency(&clockFrequency);

    auto sp = ParametersFactory::GenerateDefault();
    auto es = new CPUElasticScattering();

    static v2      tau_bounds = { 1e-13, 1e-10 };
    static cl_int2 count_bounds = { 1, 50000 };
    static v2      radius_bounds = { 1e-9, 1e-6 };
    static v2      region_bounds = { 1e-8,  1e-4 };
    static v2      extends_bounds = { 1e-7,  1e-4 };

    static v2      particle_speed_bounds = { 1e6, 1e9 };
    static v2      phi_bounds = { 0, PI2 };
    static v2      magnetic_field_bounds = { 0, 80 };
    static v2      alpha_bounds = { 0, PI/4.0 };

    static bool sync_immediate = true;
    static bool is_electron = (sp.is_clockwise == 1);
    static bool is_diag_regions = (sp.is_diag_regions == 1);
    static bool is_incoherent = (sp.is_incoherent == 1);

    const int RESULT_HISTORY_SIZE = 40;

    std::vector<float> last_results;
    last_results.resize(RESULT_HISTORY_SIZE);
    unsigned int history_index = 0;

    double last_result = 0;
    double last_result_time = 0;

    static int imp_seed = sp.impurity_seed;

    es->Compute(sp, last_result);

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {
            if (!ImGui::Begin("Elastic scattering")) {
                ImGui::End();
                continue;
            }

            static int m = sp.mode;

            ImGui::Text("Mode:"); ImGui::SameLine();
            ImGui::RadioButton("LT", &m, (int)MODE_DIR_LIFETIME); ImGui::SameLine();
            ImGui::RadioButton("PHI", &m, (int)MODE_PHI_LIFETIME); ImGui::SameLine();
            ImGui::RadioButton("S-XX", &m, (int)MODE_SIGMA_XX); ImGui::SameLine();
            ImGui::RadioButton("S-XY", &m, (int)MODE_SIGMA_XY);
            sp.mode = m;
            
            ImGui::Text("Size:"); ImGui::SameLine();
            ImGui::RadioButton("64x64", &sp.dim, 64); ImGui::SameLine();
            ImGui::RadioButton("128x128", &sp.dim, 128); ImGui::SameLine();
            ImGui::RadioButton("256x256", &sp.dim, 256); ImGui::SameLine();
            ImGui::RadioButton("1024x1024", &sp.dim, 1024);

            if (sp.mode != MODE_DIR_LIFETIME) {
                ImGui::Text("Steps:"); ImGui::SameLine();
                ImGui::RadioButton("7", &sp.integrand_steps, 7); ImGui::SameLine();
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
            ImGui::SliderScalar("Alpha", ImGuiDataType_Double, &sp.alpha, &alpha_bounds.x, &alpha_bounds.y, "%.2f");

            ImGui::Checkbox("Clockwise", &is_electron);
            sp.is_clockwise = is_electron ? 1 : 0;

            if (sp.mode != MODE_DIR_LIFETIME) {
                ImGui::Checkbox("Diag. regions", &is_diag_regions);
                ImGui::Checkbox("Incoherent", &is_incoherent);

                sp.is_diag_regions = is_diag_regions ? 1 : 0;
                sp.is_incoherent = is_incoherent ? 1 : 0;
            }

            ImGui::SliderScalar("Particle Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");
            
            ImGui::Dummy(ImVec2(0.0f, 10.0f));

            if (ImGui::CollapsingHeader("Impurity settings")) {
                ImGui::SliderInt("Count", &sp.impurity_count, count_bounds.x, count_bounds.y);
                ImGui::SliderScalar("Region", ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
                ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
                ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");
            }
            
            bool impurities_updated = ImGui::Button("New seed");
            if (impurities_updated) {
                sp.impurity_seed += 1;
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            bool update = ImGui::Button("Compute"); ImGui::SameLine();
            ImGui::Checkbox("Sync immediately", &sync_immediate);

            if (sync_immediate || update || impurities_updated) {
                QueryPerformanceCounter(&beginClock);

                double result;
                if (es->Compute(sp, result)) {
                    last_results[history_index] = result;
                    last_result = result;

                    QueryPerformanceCounter(&endClock);
                    last_result_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

                    history_index++;
                    if (history_index >= RESULT_HISTORY_SIZE)
                        history_index = 0;
                }
            }

            ImGui::Text("Mean free path: %.3e (%.2fms)", last_result, last_result_time*1000.0);
            ImGui::PlotLines("History", last_results.data(), 40, 0, 0, 0, 1e-6, ImVec2(300, 120));

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