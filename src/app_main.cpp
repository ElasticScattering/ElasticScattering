#pragma once

#include <windows.h>

#include "src/scattering/ElasticScattering.h"
#include "app_main.h"
#include "src/utils/ParametersFactory.h"
#include "src/SimulationResult.h"
#include "views/View.h"

#include <thread>
#include <chrono>

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include "implot/implot.h"


GLFWwindow* window;

ParameterView parameter_view;
TexturesView textures_view;
LogView log_view;

LARGE_INTEGER beginClock, endClock, clockFrequency;

double last_result = 0;
double last_result_time = 0;


int app_main(const InitParameters& init)
{
    SetupContexts();

    auto sp = ParametersFactory::GenerateMinimal();
    ElasticScatteringCPU es;

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput();

        {
            glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
        }

        MainLoop(es, sp);

        {
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        // Lock to 60fps.
        long time_remaining = last_result_time * 1000 - 16;
        if (time_remaining < 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(abs(time_remaining)));
        }
    }

    ShutDown();
    return 0;
}

void MainLoop(ElasticScattering &es, ScatteringParameters& sp)
{
    static bool dock_space_open = true;
    static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;

    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->GetWorkPos());
    ImGui::SetNextWindowSize(viewport->GetWorkSize());
    ImGui::SetNextWindowViewport(viewport->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
    window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;

    if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
        window_flags |= ImGuiWindowFlags_NoBackground;

    ImGui::Begin("DockSpace", &dock_space_open, window_flags);

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

        ImGui::EndMenuBar();
    }
    ImGui::PopStyleVar();

    ImGui::PopStyleVar(2);


    if (parameter_view.ShowView(sp)) {
        QueryPerformanceCounter(&beginClock);

        es.GenerateTextures(sp);
        QueryPerformanceCounter(&endClock);
        last_result_time = double(endClock.QuadPart - beginClock.QuadPart) / clockFrequency.QuadPart;
    }

    ImGui::End(); // Dockspace
}

void ProcessInput()
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void SetupContexts()
{
    if (!glfwInit()) {
        printf("GLFW init error");
        exit(-1);
    }

    glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);

    int width = 1600;
    int height = 1000;
    window = glfwCreateWindow(width, height, "Elastic Scattering", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        printf("GLFW createWindow error");
        exit(-1);
    }
    glfwMakeContextCurrent(window);

    RECT rect;
    GetClientRect(GetDesktopWindow(), &rect);
    //glfwSetWindowPos(window, rect.right / 2 - width / 2, rect.bottom / 2 - height / 2);

    glViewport(0, 0, rect.right - rect.left, rect.top - rect.bottom);

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
}

void ShutDown()
{
    ImPlot::DestroyContext();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

}

