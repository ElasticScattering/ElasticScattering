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

double last_result = 0;
std::vector<float> last_results;
const int RESULT_HISTORY_SIZE = 40;

void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void ImGuiRender(ElasticScattering &es) {

    static v2      tau_bounds = { 1e-13, 1e-10 };
    static cl_int2 count_bounds = { 1, 50000 };
    static v2      radius_bounds = { 1e-9, 1e-6 };
    static v2      region_bounds = { 1e-8,  1e-4 };
    static v2      extends_bounds = { 1e-7,  1e-4 };
    static v2      density_bounds = { 1e-7,  1e-4 };
    
    static v2      particle_speed_bounds = { 1e6, 1e9 };
    static v2      phi_bounds = { 0, PI2 };
    static v2      magnetic_field_bounds = { 0, 80 };
    static v2      alpha_bounds = { 0, PI / 4.0 };
    static v2      temperature_bounds = { 1, 300 };
    
    static bool sync_immediate = true;
    static bool is_electron = (sp.is_clockwise == 1);
    static bool is_diag_regions = (sp.is_diag_regions == 1);
    static bool is_incoherent = (sp.is_incoherent == 1);


    unsigned int history_index = 0;

    double last_result_time = 0;

    int imp_seed = sp.impurity_seed;

    //Setup dockspace
    static bool dock_space_open = true;
    static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;
    
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));

    // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
    // because it would be confusing to have two docking targets within each others.
    ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;

    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->GetWorkPos());
    ImGui::SetNextWindowSize(viewport->GetWorkSize());
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;

    // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
    // and handle the pass-thru hole, so we ask Begin() to not render a background.
    if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
        window_flags |= ImGuiWindowFlags_NoBackground;

    ImGui::Begin("DockSpace Demo", &dock_space_open, window_flags);

    ImGuiIO& io = ImGui::GetIO();
    if (io.ConfigFlags & ImGuiConfigFlags_DockingEnable)
    {
        ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
        ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
    }

    // Menu bar.
    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("Options"))
        {
            //ImGui::MenuItem("Debug information", NULL, &opt_fullscreen);
            //ImGui::MenuItem("Hide windows", NULL, &opt_fullscreen);

            ImGui::EndMenu();
        }

        ImGui::SameLine(ImGui::GetWindowWidth() - 300.0f);
        ImGui::Text("Mean free path: %.3e km? | %.2f ms", last_result, last_result_time * 1000.0);
        ImGui::EndMenuBar();
    }
    ImGui::PopStyleVar(3);

    // Elastic Scattering UI.
    bool impurities_updated = false;
    bool force_compute = false;
    bool run_iteration = false;
    {
        if (ImGui::Begin("Elastic scattering")) {
            static int m = sp.mode;
            
            {
                const char* items[] = { "One Shot", "Phi Integrated", "Sigma XX", "Sigma XY" };
                static int item_current_idx = 0;                    // Here our selection data is an index.
                const char* combo_label = items[item_current_idx];  // Label to preview before opening the combo (technically could be anything)(
                if (ImGui::BeginCombo("Mode", combo_label))
                {
                    for (int n = 0; n < IM_ARRAYSIZE(items); n++)
                    {
                        const bool is_selected = (item_current_idx == n);
                        if (ImGui::Selectable(items[n], is_selected))
                            item_current_idx = n;

                        // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
                        if (is_selected)
                            ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }

                sp.mode = item_current_idx;
            }

            {
                const char* items[] = { "64 x 64", "128 x 128", "256 x 256", "1024 x 1024" };
                const int values[]  = {  64,        128,         256,         1024 };

                static int item_current_idx = 1;
                const char* combo_label = items[item_current_idx];
                if (ImGui::BeginCombo("Size", combo_label)) {
                    for (int n = 0; n < IM_ARRAYSIZE(items); n++) {
                        const bool is_selected = (item_current_idx == n);
                        if (ImGui::Selectable(items[n], is_selected))
                            item_current_idx = n;

                        if (is_selected)
                            ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }

                sp.dim = values[item_current_idx];
            }

            if (sp.mode != MODE_DIR_LIFETIME) {
                const char* items[] = { "7", "13", "25", "49", "99" };
                const int values[]  = {  7,   13,   25,   49,   99 };

                static int item_current_idx = 1;
                const char* combo_label = items[item_current_idx];
                if (ImGui::BeginCombo("Steps", combo_label)) {
                    for (int n = 0; n < IM_ARRAYSIZE(items); n++) {
                        const bool is_selected = (item_current_idx == n);
                        if (ImGui::Selectable(items[n], is_selected))
                            item_current_idx = n;

                        if (is_selected)
                            ImGui::SetItemDefaultFocus();
                    }
                    ImGui::EndCombo();
                }

                sp.integrand_steps = values[item_current_idx];
            } 
            
            ImGui::Checkbox("Clockwise", &is_electron);
            sp.is_clockwise = is_electron ? 1 : 0;

            if (sp.mode != MODE_DIR_LIFETIME) {
                ImGui::Checkbox("Diagonal regions", &is_diag_regions);
                ImGui::Checkbox("Incoherent", &is_incoherent);

                sp.is_diag_regions = is_diag_regions ? 1 : 0;
                sp.is_incoherent = is_incoherent ? 1 : 0;
            }
            ImGui::SliderScalar("Alpha", ImGuiDataType_Double, &sp.alpha, &alpha_bounds.x, &alpha_bounds.y, "%.2f");

            ImGui::Dummy(ImVec2(0.0f, 15.0f));

            if (is_incoherent) ImGui::SliderScalar("Temperature (K)", ImGuiDataType_Double, &sp.temperature, &temperature_bounds.x, &temperature_bounds.y, "%.1f");
            else               ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%.2e");
            
            ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%.2f");
            if (sp.mode == MODE_DIR_LIFETIME)
                ImGui::SliderScalar("Phi", ImGuiDataType_Double, &sp.phi, &phi_bounds.x, &phi_bounds.y, "%.2f", 1.0f);

            ImGui::SliderScalar("Particle Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");

            ImGui::Dummy(ImVec2(0.0f, 15.0f));

            ImGui::Text("Impurity settings");
            ImGui::SliderInt("Count", &sp.impurity_count, count_bounds.x, count_bounds.y);
            ImGui::SliderScalar("Region",  ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
            ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
            ImGui::SliderScalar("Density", ImGuiDataType_Double, &sp.impurity_density, &density_bounds.x, &density_bounds.y, "%.2e");
            ImGui::SliderScalar("Radius",  ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");

            impurities_updated = ImGui::Button("New seed");
            if (impurities_updated) {
                sp.impurity_seed += 1;
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            force_compute = ImGui::Button("Compute"); ImGui::SameLine();
            ImGui::Checkbox("Sync immediately", &sync_immediate);

            run_iteration = ImGui::Button("Iterate over B: 0 .. 50");
        }
        ImGui::End();
    }

    // Disable everything and run this iteration, draw frame until 50 frames have passed.
    if (run_iteration) {
        history_index = 0;
        double result, result2, result3, result4;

        for (int i = 0; i < 50; i++) {
            sp.magnetic_field = i;
            
            sp.is_incoherent = false;
            sp.mode = MODE_SIGMA_XX;
            es.Compute(sp, result);

            sp.mode = MODE_SIGMA_XX;
            es.Compute(sp, result);

            sp.is_incoherent = true;
            sp.mode = MODE_SIGMA_XY;
            es.Compute(sp, result);

            sp.mode = MODE_SIGMA_XY;
            es.Compute(sp, result);


            last_results[history_index] = result;

            history_index++;
            if (history_index >= RESULT_HISTORY_SIZE)
                history_index = 0;
        }

    } else if (sync_immediate || force_compute || impurities_updated) {
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

