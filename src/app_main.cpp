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

LARGE_INTEGER beginClock, endClock, clockFrequency;

SimulationParameters sp;

v2      tau_bounds = { 1e-13, 1e-10 };
cl_int2 count_bounds = { 1, 50000 };
v2      radius_bounds = { 1e-9, 1e-6 };
v2      region_bounds = { 1e-8,  1e-4 };
v2      extends_bounds = { 1e-7,  1e-4 };
v2      density_bounds = { 1e-7,  1e-4 };

v2      particle_speed_bounds = { 1e6, 1e9 };
v2      phi_bounds = { 0, PI2 };
v2      magnetic_field_bounds = { 0, 80 };
v2      alpha_bounds = { 0, PI / 4.0 };
v2      temperature_bounds = { 1, 300 };

bool sync_immediate = true;
bool is_electron = (sp.is_clockwise == 1);
bool is_diag_regions = (sp.is_diag_regions == 1);
bool is_incoherent = (sp.is_incoherent == 1);

const int RESULT_HISTORY_SIZE = 40;

std::vector<float> last_results;
unsigned int history_index = 0;

double last_result = 0;
double last_result_time = 0;

int imp_seed = sp.impurity_seed;


void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void ImGuiRender(ElasticScattering &es) {

    //Setup dockspace
    static bool opt_fullscreen = true;
    static bool dock_space_open = true;
    static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;
    
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));

    // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
    // because it would be confusing to have two docking targets within each others.
    ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
    if (opt_fullscreen)
    {
        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->GetWorkPos());
        ImGui::SetNextWindowSize(viewport->GetWorkSize());
        ImGui::SetNextWindowViewport(viewport->ID);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
        window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
        window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
    }
    else
    {
        dockspace_flags &= ~ImGuiDockNodeFlags_PassthruCentralNode;
    }

    // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
    // and handle the pass-thru hole, so we ask Begin() to not render a background.
    if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
        window_flags |= ImGuiWindowFlags_NoBackground;

    ImGui::Begin("DockSpace Demo", &dock_space_open, window_flags);

    if (opt_fullscreen)
        ImGui::PopStyleVar(2);

    ImGuiIO& io = ImGui::GetIO();
    if (io.ConfigFlags & ImGuiConfigFlags_DockingEnable)
    {
        ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
        ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
    }

    // Menu bar.
    if (ImGui::BeginMainMenuBar())
    {
        if (ImGui::BeginMenu("Options"))
        {
            ImGui::MenuItem("Fullscreen", NULL, &opt_fullscreen);
            //ImGui::MenuItem("Debug information", NULL, &opt_fullscreen);
            //ImGui::MenuItem("Hide windows", NULL, &opt_fullscreen);

            ImGui::EndMenu();
        }

        ImGui::SameLine(ImGui::GetWindowWidth() - 320.0f);
        ImGui::Text("Mean free path: %.3e km? | %.2f ms", last_result, last_result_time * 1000.0);
        ImGui::EndMainMenuBar();
    }
    ImGui::PopStyleVar();
  

    // Elastic Scattering UI.
    bool impurities_updated = false;
    bool force_compute = false;
    {
        if (ImGui::Begin("Elastic scattering")) {
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

            ImGui::SliderScalar("Temperature (K)", ImGuiDataType_Double, &sp.temperature, &temperature_bounds.x, &temperature_bounds.y, "%.2e");
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
                ImGui::SliderScalar("Density", ImGuiDataType_Double, &sp.impurity_density, &density_bounds.x, &density_bounds.y, "%.2e");

                ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");
            }

            impurities_updated = ImGui::Button("New seed");
            if (impurities_updated) {
                sp.impurity_seed += 1;
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            force_compute = ImGui::Button("Compute"); ImGui::SameLine();
            ImGui::Checkbox("Sync immediately", &sync_immediate);
        }
        ImGui::End();
    }

    if (sync_immediate || force_compute || impurities_updated) {
        QueryPerformanceCounter(&beginClock);

        double result;
        if (es.Compute(sp, result)) {
            last_results[history_index] = result;
            last_result = result;

            QueryPerformanceCounter(&endClock);
            last_result_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

            history_index++;
            if (history_index >= RESULT_HISTORY_SIZE)
                history_index = 0;
        }
    }
    
    //ImGui::PlotLines("History", last_results.data(), 40, 0, 0, 0, 1e-6, ImVec2(300, 120));


    // Output texture view
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("Output");

    auto texture = es.GetTextureID();

    auto viewportSize = ImGui::GetContentRegionAvail();
    float min_dim = min(viewportSize.x, viewportSize.y);
    ImGui::Image((void*)texture, ImVec2(min_dim,min_dim));

    ImGui::PopStyleVar();
    ImGui::End();

    ImGui::End(); // Dockspace
}

int app_main(int argc, char** argv)
{
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);

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
    //glfwSetWindowPos(window, rect.right / 2 - width / 2, rect.bottom / 2 - height / 2);

    glViewport(0, 0, rect.right-rect.left, rect.top - rect.bottom);

    GLenum error = glewInit();
    if (GLEW_OK != error)
    {
        printf("ERROR: %s", glewGetErrorString(error));
        exit(EXIT_FAILURE);
    }

    ImGui::CreateContext();
    ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");

    QueryPerformanceFrequency(&clockFrequency);

    sp = ParametersFactory::GenerateDefault();
    CPUElasticScattering es;

    last_results.resize(RESULT_HISTORY_SIZE);

    es.Compute(sp, last_result);

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGuiRender(es);

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

