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

#include "implot/implot.h"

LARGE_INTEGER beginClock, endClock, clockFrequency;

SimulationParameters sp;

double last_result = 0;
double last_result_time = 0;

const int RESULT_HISTORY_SIZE = 50;
std::vector<double> xs;
std::vector<double> xs_rho;

std::vector<double> last_iteration_results_xx_coherent;
std::vector<double> last_iteration_results_xy_coherent;
std::vector<double> last_iteration_results_xx_incoherent;
std::vector<double> last_iteration_results_xy_incoherent;

std::vector<double> last_iteration_results_rho_coherent;
std::vector<double> last_iteration_results_rho_incoherent;

std::vector<double> last_iteration_results_rho;

std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax;
std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax2;
std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax3;
std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax4;

std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax_rho;

bool first_iteration = true;

void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void ImGuiRender(ElasticScattering &es) {

    static v2      tau_bounds = { 1e-13, 1e-10 };
    static v2      temperature_bounds = { 0.001, 0.99 };
    
    static v2      radius_bounds  = { 5e-8, 2e-6 };
    static v2      region_bounds  = { 1e-6,  5e-4 };
    static v2      extends_bounds = { 1e-5,  5e-4 };
    static v2      density_bounds = { 1e9,  1e12 };
    
    static v2      particle_speed_bounds = { 1e6, 1e9 };
    static v2      phi_bounds = { 0, PI2 };
    static v2      magnetic_field_bounds = { 0, 80 };
    static v2      alpha_bounds = { 0, PI / 4.0 };
    
    static bool is_electron = (sp.is_clockwise == 1);
    static bool is_diag_regions = (sp.is_diag_regions == 1);
    static bool is_incoherent = (sp.is_incoherent == 1);
    
    static bool sync_immediate = true;
    static bool interactive = true;

    unsigned int history_index = 0;

    int imp_seed = sp.impurity_seed;

    //Setup dockspace
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

    ImGui::PopStyleVar(2);

    // Elastic Scattering UI.
    bool impurities_updated = false;
    bool force_compute = false;
    bool run_iteration = false;
    {
        if (ImGui::Begin("Elastic scattering")) {
            static int m = sp.mode;
            
            {
                const char* items[] = { "Interactive", "Graph" };
                static int item_current_idx = 0;                    // Here our selection data is an index.
                const char* combo_label = items[item_current_idx];  // Label to preview before opening the combo (technically could be anything)(
                if (ImGui::BeginCombo("Program", combo_label))
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

                interactive = (item_current_idx == 0);
            }

            if (interactive) {
                const char* items[] = { "One Path", "Phi Integrated", "Sigma XX", "Sigma XY"};
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

            if (!interactive || sp.mode != MODE_DIR_LIFETIME) {
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

            if (!interactive || sp.mode != MODE_DIR_LIFETIME) {
                ImGui::Checkbox("Diagonal regions", &is_diag_regions);
                ImGui::Checkbox("Incoherent", &is_incoherent);

                sp.is_diag_regions = is_diag_regions ? 1 : 0;
                sp.is_incoherent = is_incoherent ? 1 : 0;
            }

            ImGui::SliderScalar("Alpha", ImGuiDataType_Double, &sp.alpha, &alpha_bounds.x, &alpha_bounds.y, "%.2f");

            ImGui::Dummy(ImVec2(0.0f, 15.0f));

            if (is_incoherent) ImGui::SliderScalar("Temperature (K)", ImGuiDataType_Double, &sp.temperature, &temperature_bounds.x, &temperature_bounds.y, "%.3f");
            else               ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%.2e");
            
            ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%.2f");
            if (interactive && sp.mode == MODE_DIR_LIFETIME)
                ImGui::SliderScalar("Phi", ImGuiDataType_Double, &sp.phi, &phi_bounds.x, &phi_bounds.y, "%.2f", 1.0f);

            ImGui::SliderScalar("Particle Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");

            ImGui::Dummy(ImVec2(0.0f, 15.0f));
            
            ImGui::Text("Impurity settings");
            ImGui::SliderScalar("Region",  ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
            ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
            ImGui::SliderScalar("Density", ImGuiDataType_Double, &sp.impurity_density, &density_bounds.x, &density_bounds.y, "%.2e");
            ImGui::SliderScalar("Radius",  ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");

            impurities_updated = ImGui::Button("New seed");
            if (impurities_updated) {
                sp.impurity_seed += 1;
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            if (interactive) {
                force_compute = ImGui::Button("Compute"); ImGui::SameLine();
                ImGui::Checkbox("Sync immediately", &sync_immediate);
            }
            else {
                run_iteration = ImGui::Button("Generate graphs");
            }
        }
        ImGui::End();
    }

    if (interactive && (sync_immediate || force_compute || impurities_updated)) {
        QueryPerformanceCounter(&beginClock);

        double result;
        if (es.Compute(sp, result)) {
            last_result = result;

            QueryPerformanceCounter(&endClock);
            last_result_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
        }
    }
    else if (run_iteration) { // Graph / Not interactive
        run_iteration = false;
        first_iteration = false;
        history_index = 0;
        double result, result2, result3, result4, result_rho, result_rho_;

        double last_radius = sp.impurity_radius;
        double last_count  = sp.impurity_count;
        double last_field  = sp.magnetic_field;
        int    last_coherent = sp.is_incoherent;
        int    last_mode = sp.mode;

        sp.impurity_radius = 1e-18;
        sp.impurity_count = 1;

        for (int i = 0; i < RESULT_HISTORY_SIZE; i++) {
            sp.magnetic_field = i;

            sp.is_incoherent = false;
            sp.mode = MODE_SIGMA_XX;
            es.Compute(sp, result);

            sp.mode = MODE_SIGMA_XY;
            es.Compute(sp, result2);

            sp.is_incoherent = true;
            sp.mode = MODE_SIGMA_XX;
            es.Compute(sp, result3);

            sp.mode = MODE_SIGMA_XY;
            es.Compute(sp, result4);

            printf("%i %.12e %.12e %.12e %.12e\n", i, result, result2, result3, result4);
            last_iteration_results_xx_coherent[i] = result / 1e8;
            last_iteration_results_xy_coherent[i] = result2 / 1e8;
            last_iteration_results_xx_incoherent[i] = result3 / 1e8;
            last_iteration_results_xy_incoherent[i] = result4 / 1e8;

            double sxx = last_iteration_results_xx_incoherent[i] + last_iteration_results_xx_coherent[i];
            double sxy = last_iteration_results_xy_incoherent[i] + last_iteration_results_xy_coherent[i];
            last_iteration_results_rho_incoherent[i] = sxx / ((sxx * sxx) + (sxy * sxy)); //incoherent
        }

        for (int i = 1; i < RESULT_HISTORY_SIZE; i++) {
            double Y = last_iteration_results_rho_incoherent[i] - last_iteration_results_rho_incoherent[i - 1]; //stappen van 1

            last_iteration_results_rho[i-1] = Y;
        }

        printf("\n");

        minmax  = std::minmax_element(last_iteration_results_xx_coherent.begin(), last_iteration_results_xx_coherent.end());
        minmax2 = std::minmax_element(last_iteration_results_xy_coherent.begin(), last_iteration_results_xy_coherent.end());
        minmax3 = std::minmax_element(last_iteration_results_xx_incoherent.begin(), last_iteration_results_xx_incoherent.end());
        minmax4 = std::minmax_element(last_iteration_results_xy_incoherent.begin(), last_iteration_results_xy_incoherent.end());

        minmax_rho = std::minmax_element(last_iteration_results_rho.begin(), last_iteration_results_rho.end());

        sp.impurity_radius = last_radius;
        sp.impurity_count  = last_count;
        sp.magnetic_field = last_field;
        sp.is_incoherent = last_coherent;
        sp.mode = last_mode;
    }

    //ImPlot::ShowDemoWindow(0);
    
    ImGui::Begin("Iteration");
    if (!first_iteration) {
        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax3.first, *minmax3.second);
        if (ImPlot::BeginPlot("Sigma XX Incoherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XX Incoherent", xs.data(), last_iteration_results_xx_incoherent.data(), RESULT_HISTORY_SIZE);
            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax4.first, *minmax4.second);
        if (ImPlot::BeginPlot("Sigma XY Incoherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XY Incoherent", xs.data(), last_iteration_results_xy_incoherent.data(), RESULT_HISTORY_SIZE);

            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax.first, *minmax.second);
        if (ImPlot::BeginPlot("Sigma XX Coherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XX Coherent", xs.data(), last_iteration_results_xx_coherent.data(), RESULT_HISTORY_SIZE);
            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax2.first, *minmax2.second);
        if (ImPlot::BeginPlot("Sigma XY Coherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XY Coherent", xs.data(), last_iteration_results_xy_coherent.data(), RESULT_HISTORY_SIZE);

            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, (RESULT_HISTORY_SIZE-1), *minmax_rho.first, *minmax_rho.second);
        if (ImPlot::BeginPlot("Rho", "Magnetic field (tesla)", "Afgeleide rho")) {
            ImPlot::PlotLine("Rho", xs_rho.data(), last_iteration_results_rho.data(), RESULT_HISTORY_SIZE-1);

            ImPlot::EndPlot();
        }
    }
    ImGui::End();

    // Output texture view
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("Output");

    auto texture = es.GetTextureID();

    auto viewportSize = ImGui::GetContentRegionAvail();
    float min_dim = min(viewportSize.x, viewportSize.y);
    ImGui::Image((void*)texture, ImVec2(min_dim,min_dim));



    ImGui::PopStyleVar();
    ImGui::End();

    // Menu bar.
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));

    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("Options"))
        {
            //ImGui::MenuItem("Debug information", NULL, &opt_fullscreen);
            //ImGui::MenuItem("Hide windows", NULL, &opt_fullscreen);

            ImGui::EndMenu();
        }

        ImGui::SameLine(ImGui::GetWindowWidth() - 300.0f);
        ImGui::Text("Mean free path: %.3e m | (%.3f ms)", last_result, (last_result_time * 1000.0));
        ImGui::EndMenuBar();
    }
    ImGui::PopStyleVar();


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

    ImPlot::CreateContext();

    QueryPerformanceFrequency(&clockFrequency);

    sp = ParametersFactory::GenerateDefault();
    GPUElasticScattering es;

    last_iteration_results_xx_coherent.resize(RESULT_HISTORY_SIZE);
    last_iteration_results_xy_coherent.resize(RESULT_HISTORY_SIZE);
    last_iteration_results_xx_incoherent.resize(RESULT_HISTORY_SIZE);
    last_iteration_results_xy_incoherent.resize(RESULT_HISTORY_SIZE);
    last_iteration_results_rho_incoherent.resize(RESULT_HISTORY_SIZE);
    last_iteration_results_rho.resize(RESULT_HISTORY_SIZE-1);

    xs.resize(RESULT_HISTORY_SIZE);
    for (int i = 0; i < RESULT_HISTORY_SIZE; i++) {
        xs[i] = i;
    }

    xs_rho.resize(RESULT_HISTORY_SIZE-1);
    for (int i = 0; i < RESULT_HISTORY_SIZE-1; i++) {
        xs_rho[i] = i;
    }

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

    ImPlot::DestroyContext();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

