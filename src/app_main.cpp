#pragma once
#ifndef NO_WINDOW
#include <windows.h>
/*
#include "src/scattering/ElasticScattering.h"
#include "app_main.h"
#include "utils/OpenCLUtils.h"
#include "utils/ErrorMacros.h"
#include "src/ParametersFactory.h"
#include "src/escl/constants.h"
#include "src/SimulationResult.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include "math.h"
#include <thread>
#include <chrono>

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "implot/implot.h"

LARGE_INTEGER beginClock, endClock, clockFrequency;

ScatteringParameters sp;

double last_result = 0;
double last_result_time = 0;

bool first_iteration = true;
bool sync_immediate = true;

v2 tau_bounds            = { 1e-13, 1e-10 };
v2 temperature_bounds    = { 0, 10 };

v2 radius_bounds         = { 5e-8, 2e-6 };
v2 region_bounds         = { 1e-6,  5e-4 };
v2 extends_bounds        = { 1e-5,  5e-4 };
v2 density_bounds        = { 1e9,  1e13 };

v2 particle_speed_bounds = { 1e6, 1e9 };
v2 phi_bounds            = { 0, PI2 };
v2 magnetic_field_bounds = { 0, 80 };
v2 alpha_bounds          = { 0, PI / 4.0 };

typedef struct {
    bool compute_requested;
} Action;

void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

Action ShowWindowParameters() {
    bool is_electron = (sp.is_clockwise == 1);
    bool is_diag_regions = (sp.is_diag_regions == 1);
    bool is_incoherent = (sp.is_incoherent == 1);

    Action action;
    action.compute_requested = false;

    if (ImGui::Begin("Elastic scattering")) {
        static int m = sp.mode;

        {
            const char* items[] = { "One Path", "Phi Integrated", "Sigma XX", "Sigma XY" };
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
            const int values[] = { 64,        128,         256,         1024 };

            static int item_current_idx = 0;
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
            const int values[] = { 7,   13,   25,   49,   99 };

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

        ImGui::Checkbox("Diagonal regions", &is_diag_regions);
        sp.is_diag_regions = is_diag_regions ? 1 : 0;

        if (sp.mode == MODE_DIR_LIFETIME) {
            ImGui::Checkbox("Incoherent", &is_incoherent);

            sp.is_incoherent = is_incoherent ? 1 : 0;
        }

        ImGui::SliderScalar("Alpha", ImGuiDataType_Double, &sp.alpha, &alpha_bounds.x, &alpha_bounds.y, "%.2f");

        ImGui::Dummy(ImVec2(0.0f, 15.0f));

        if (is_incoherent) ImGui::SliderScalar("Temperature (K)", ImGuiDataType_Double, &sp.temperature, &temperature_bounds.x, &temperature_bounds.y, "%.3f");
        else               ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%.2e");

        ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%.2f");
            
        if (sp.mode == MODE_DIR_LIFETIME)
            ImGui::SliderScalar("Phi", ImGuiDataType_Double, &sp.phi, &phi_bounds.x, &phi_bounds.y, "%.2f", 1.0f);

        ImGui::SliderScalar("Particle Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");

        ImGui::Dummy(ImVec2(0.0f, 15.0f));

        ImGui::Text("Impurity settings");
        ImGui::SliderScalar("Region", ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
        ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
        ImGui::SliderScalar("Density", ImGuiDataType_Double, &sp.impurity_density, &density_bounds.x, &density_bounds.y, "%.2e");
        ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");

        if (ImGui::Button("New seed")) {
            sp.impurity_seed += 1;
            action.compute_requested = true;
        }
        
        ImGui::Dummy(ImVec2(0.0f, 20.0f));

        action.compute_requested |= ImGui::Button("Compute"); ImGui::SameLine();
        ImGui::Checkbox("Sync immediately", &sync_immediate);
    }
    ImGui::End();

    return action;
}

void ImGuiRender(ElasticScattering &es) {
    //Setup dockspace
    {
        static bool dock_space_open = true;
        static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;

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
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
        if (ImGui::BeginMenuBar())
        {
            ImGui::SameLine(ImGui::GetWindowWidth() - 300.0f);
            ImGui::Text("Mean free path: %.3e m | (%.3f ms)", last_result, (last_result_time * 1000.0));
            //printf("Mean free path: %.3e m | (%.3f ms)\n", last_result, (last_result_time * 1000.0));

            ImGui::EndMenuBar();
        }
        ImGui::PopStyleVar();

        ImGui::PopStyleVar(2);
    }


    Action action = ShowWindowParameters();

    if (sync_immediate || action.compute_requested) {
        QueryPerformanceCounter(&beginClock);

        if (es.Compute(sp, last_result)) {
            QueryPerformanceCounter(&endClock);
            last_result_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
            //std::cout << last_result_time << " " << sp.impurity_count << std::endl;
        }
    }

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

int app_main(const InitParameters& init)
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

    ImPlot::CreateContext();

    QueryPerformanceFrequency(&clockFrequency);

    sp = ParametersFactory::GenerateMinimal();
    GPUElasticScattering es(init);

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

        long time_remaining = last_result_time * 1000 - 16;
        if (time_remaining < 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(abs(time_remaining)));
        }
    }

    ImPlot::DestroyContext();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
*/
#endif // NO_WINDOW