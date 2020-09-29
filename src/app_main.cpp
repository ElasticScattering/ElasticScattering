#pragma once
#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "ElasticScattering.h"
#include "app_main.h"
#include "utils/OpenCLUtils.h"
#include "utils/ErrorMacros.h"
#include "src/ParametersFactory.h"
#include "src/Logger.h"
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

SimulationParameters sp;

double last_result = 0;
double last_result_time = 0;

bool first_iteration = true;
bool sync_immediate = true;

Logger logger;

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
    bool interactive;
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
    action.interactive = false;
    action.compute_requested = false;

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

            action.interactive = (item_current_idx == 0);
        }

        if (action.interactive) {
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

        if (!action.interactive || sp.mode != MODE_DIR_LIFETIME) {
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

        if (action.interactive && sp.mode == MODE_DIR_LIFETIME) {
            ImGui::Checkbox("Incoherent", &is_incoherent);

            sp.is_incoherent = is_incoherent ? 1 : 0;
        }

        ImGui::SliderScalar("Alpha", ImGuiDataType_Double, &sp.alpha, &alpha_bounds.x, &alpha_bounds.y, "%.2f");

        ImGui::Dummy(ImVec2(0.0f, 15.0f));

        if (is_incoherent) ImGui::SliderScalar("Temperature (K)", ImGuiDataType_Double, &sp.temperature, &temperature_bounds.x, &temperature_bounds.y, "%.3f");
        else               ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%.2e");

        if (action.interactive) {
            ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%.2f");
            
            if (sp.mode == MODE_DIR_LIFETIME)
                ImGui::SliderScalar("Phi", ImGuiDataType_Double, &sp.phi, &phi_bounds.x, &phi_bounds.y, "%.2f", 1.0f);
        }

        ImGui::SliderScalar("Particle Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");

        ImGui::Dummy(ImVec2(0.0f, 15.0f));

        ImGui::Text("Impurity settings");
        ImGui::SliderScalar("Region", ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
        ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
        ImGui::SliderScalar("Density", ImGuiDataType_Double, &sp.impurity_density, &density_bounds.x, &density_bounds.y, "%.2e");
        ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");

        if (action.interactive) {
            if (ImGui::Button("New seed")) {
                sp.impurity_seed += 1;
                action.compute_requested = true;
            }
        }
        
        ImGui::Dummy(ImVec2(0.0f, 20.0f));

        if (action.interactive) {
            action.compute_requested |= ImGui::Button("Compute"); ImGui::SameLine();
            ImGui::Checkbox("Sync immediately", &sync_immediate);
        }
        else {
            action.compute_requested |= ImGui::Button("Generate graphs (magnetic field)");
        }
    }
    ImGui::End();

    return action;
}

void ComputeSimulation(ElasticScattering& es, SimulationResult& sr)
{
    first_iteration = false;


    /*
    std::vector<double> zs{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                    1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
    sr.n_runs = zs.size();
    */

    sr.n_runs = 20;
    sr.xs.resize(sr.n_runs);
    sr.xs_temperature.resize(sr.n_runs);

    sr.results_xx.resize(sr.n_runs);
    sr.results_xxi.resize(sr.n_runs);
    sr.delta_xxi.resize(sr.n_runs);
    sr.results_xy.resize(sr.n_runs);
    sr.results_xyi.resize(sr.n_runs);

    sr.iterations_per_run = 3;
    sr.x_is_temperature = false;

    v2 range = { 0.01, 40 };
    double step_size = (range.y - range.x) / sr.n_runs;

    double coherent_tau = sp.tau;
    bool run_incoherent = sp.alpha > 0.000001;
    bool run_coherent = abs(sp.alpha - (PI / 4)) > 0.000001;

    for (int i = 0; i < sr.n_runs; i++) {
        sp.magnetic_field = range.x + step_size * i;
        sr.xs[i] = sp.magnetic_field;
        sr.xs_temperature[i] = sp.temperature;

        // SIGMA XX
        {
            sp.mode = MODE_SIGMA_XX;

            double total_sigma_xx = 0;
            double total_sigma_xxi = 0;
            double total_sigma_xx_sq = 0;

            for (int j = 0; j < sr.iterations_per_run; j++) {
                sp.impurity_seed = 1123 + j * 831;
                double result, resulti;

                if (run_incoherent) {
                    sp.is_incoherent = 1;
                    es.Compute(sp, resulti);
                }
                else {
                    resulti = 0;
                }

                if (run_coherent) {
                    sp.is_incoherent = 0;
                    sp.tau = coherent_tau;
                    es.Compute(sp, result);
                }
                else {
                    result = 0;
                }

                double sxx = result / 1e8;
                double sxx_i = resulti / 1e8;

                total_sigma_xx += sxx;
                total_sigma_xxi += sxx_i;
                total_sigma_xx_sq += (sxx_i + sxx) * (sxx_i + sxx);
            }

            {
                sr.results_xxi[i] = (total_sigma_xxi) / (double)(sr.iterations_per_run);
                sr.results_xx[i] = (total_sigma_xx) / (double)(sr.iterations_per_run);


                double sxx_sq_exp = total_sigma_xx_sq / (double)(sr.iterations_per_run) + 1e-15;
                double sxx_exp = (total_sigma_xx + total_sigma_xxi) / (double)(sr.iterations_per_run);

                double sxx_std = sqrt((sxx_sq_exp - sxx_exp * sxx_exp) / (double)(sr.iterations_per_run-1));
                sr.delta_xxi[i] = sxx_std / sxx_exp;
            }
        }

        // SIGMA XY
        {
            sp.mode = MODE_SIGMA_XY;

            double total_sigma_xy = 0;
            double total_sigma_xyi = 0;
            double total_sigma_xy_sq = 0;

            for (int j = 0; j < sr.iterations_per_run; j++) {
                sp.impurity_seed = 1123 + j * 831;
                double result, resulti;

                if (run_incoherent) {
                    sp.is_incoherent = 1;
                    es.Compute(sp, resulti);
                }
                else {
                    resulti = 0;
                }

                if (run_coherent) {
                    sp.is_incoherent = 0;
                    sp.tau = coherent_tau;
                    es.Compute(sp, result);
                }
                else {
                    result = 0;
                }

                double sxy = result / 1e8;
                double sxy_i = resulti / 1e8;

                total_sigma_xy += sxy;
                total_sigma_xyi += sxy_i;
            }

            {
                sr.results_xy[i] = total_sigma_xy / (double)sr.iterations_per_run;
                sr.results_xyi[i] = total_sigma_xyi / (double)sr.iterations_per_run;
            }
        }

        printf("Progress %d/%d\n", i, sr.n_runs);
    }
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

    if (action.interactive && (sync_immediate || action.compute_requested)) {
        QueryPerformanceCounter(&beginClock);

        if (es.Compute(sp, last_result)) {
            QueryPerformanceCounter(&endClock);
            last_result_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
            //std::cout << last_result_time << " " << sp.impurity_count << std::endl;
        }
    }
    else if (action.compute_requested) { // Graph / Not interactive
        SimulationParameters old_sp = sp;

#if 1
        sp.impurity_density = 5.34e14;
        sp.impurity_radius = 1.11e-8;
        sp.region_extends = 1e-6;
        sp.region_size = 4e-6;
        sp.alpha = 0.3;
        sp.dim = 128;
        sp.tau = 1e-11;
#else
        sp.impurity_density = 5.34e12;
        sp.impurity_radius = 1.11e-7;
        sp.region_extends = 1e-5;
        sp.region_size = 4e-5;
        sp.alpha = PI/4;
        sp.dim = 128;
#endif

        const std::vector<double> zs{15, 60 };

        //const std::vector<double> zs{ 1 };

        for (int i = 0; i < zs.size(); i++) {
            sp.temperature = zs[i];

            QueryPerformanceCounter(&beginClock);

            SimulationResult result;
            ComputeSimulation(es, result);
            QueryPerformanceCounter(&endClock);
            result.time_elapsed = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;

            logger.LogResult(sp, result);
            printf("Simulation completed %d/%d\n", i+1, zs.size());
        }
     
        sp = old_sp;
    }

#if 0
    //std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax_xx;
    //std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax_xy;
    //std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax_xxi;
    //std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax_xyi;

    //std::pair<std::vector<double>::iterator, std::vector<double>::iterator> minmax_rho;

    auto minmax_xx = std::minmax_element(sr.results_xx.begin(), sr.results_xx.end());
    auto minmax_xy = std::minmax_element(sr.results_xy.begin(), sr.results_xy.end());
    auto minmax_xxi = std::minmax_element(sr.results_xxi.begin(), sr.results_xxi.end());
    auto minmax_xyi = std::minmax_element(sr.results_xyi.begin(), sr.results_xyi.end());

    auto minmax_rho = std::minmax_element(sr.results_rho_deriv.begin(), sr.results_rho_deriv.end());

    /*
    ImGui::Begin("Iteration");
    if (!first_iteration) {
        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax_xxi.first, *minmax_xxi.second);
        if (ImPlot::BeginPlot("Sigma XX Incoherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XX Incoherent", xs.data(), results_xxi.data(), RESULT_HISTORY_SIZE);
            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax_xyi.first, *minmax_xyi.second);
        if (ImPlot::BeginPlot("Sigma XY Incoherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XY Incoherent", xs.data(), results_xyi.data(), RESULT_HISTORY_SIZE);

            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax_xx.first, *minmax_xx.second);
        if (ImPlot::BeginPlot("Sigma XX Coherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XX Coherent", xs.data(), results_xx.data(), RESULT_HISTORY_SIZE);
            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, RESULT_HISTORY_SIZE, *minmax_xy.first, *minmax_xy.second);
        if (ImPlot::BeginPlot("Sigma XY Coherent", "Magnetic field (tesla)", "Mean free path (km)")) {
            ImPlot::PlotLine("Sigma XY Coherent", xs.data(), results_xy.data(), RESULT_HISTORY_SIZE);

            ImPlot::EndPlot();
        }

        ImPlot::SetNextPlotLimits(0, (RESULT_HISTORY_SIZE-1), *minmax_rho.first, *minmax_rho.second);
        if (ImPlot::BeginPlot("Rho", "Magnetic field (tesla)", "Afgeleide rho")) {
            ImPlot::PlotLine("Rho", xs_rho.data(), results_rho_deriv.data(), RESULT_HISTORY_SIZE-1);

            ImPlot::EndPlot();
        }
    }
    ImGui::End();
    */
#endif

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

