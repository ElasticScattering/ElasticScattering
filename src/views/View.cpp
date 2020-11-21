
#include "View.h"

#include "src/scattering/escl/constants.h"

#include <GL/glew.h>
#include <GL/glfw3.h>

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

bool sync_immediate = true;

v2 tau_bounds = { 1e-13, 1e-10 };
v2 temperature_bounds = { 0, 10 };

v2 radius_bounds = { 5e-8, 2e-6 };
v2 region_bounds = { 1e-6,  5e-4 };
v2 extends_bounds = { 1e-5,  5e-4 };
v2 density_bounds = { 1e9,  1e13 };

v2 particle_speed_bounds = { 1e6, 1e9 };
v2 phi_bounds = { 0, PI2 };
v2 magnetic_field_bounds = { 0, 80 };
v2 alpha_bounds = { 0, PI / 4.0 };

bool ParameterView::ShowView(ScatteringParameters& sp)
{
    bool is_electron = (sp.is_clockwise == 1);
    bool is_diag_regions = (sp.is_diag_regions == 1);
    bool is_incoherent = (sp.is_incoherent == 1);

    bool compute_requested = false;
    bool parameters_changed = false;

    if (ImGui::Begin("Elastic scattering")) {
        {
            const char* items[] = { "64 x 64", "128 x 128", "256 x 256", "1024 x 1024" };
            const int values[]  = {  64,        128,         256,         1024 };

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

        {
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
        ImGui::Checkbox("Diagonal regions", &is_diag_regions);
        ImGui::Checkbox("Incoherent", &is_incoherent);

        sp.is_clockwise = is_electron ? 1 : 0;
        sp.is_diag_regions = is_diag_regions ? 1 : 0;
        sp.is_incoherent = is_incoherent ? 1 : 0;

        ImGui::SliderScalar("Alpha", ImGuiDataType_Double, &sp.alpha, &alpha_bounds.x, &alpha_bounds.y, "%.2f");

        ImGui::Dummy(ImVec2(0.0f, 15.0f));

        if (is_incoherent) ImGui::SliderScalar("Temperature (K)", ImGuiDataType_Double, &sp.temperature, &temperature_bounds.x, &temperature_bounds.y, "%.3f");
        else               ImGui::SliderScalar("Tau", ImGuiDataType_Double, &sp.tau, &tau_bounds.x, &tau_bounds.y, "%.2e");

        ImGui::SliderScalar("B", ImGuiDataType_Double, &sp.magnetic_field, &magnetic_field_bounds.x, &magnetic_field_bounds.y, "%.2f");

        ImGui::SliderScalar("Particle Speed", ImGuiDataType_Double, &sp.particle_speed, &particle_speed_bounds.x, &particle_speed_bounds.y, "%.2e");

        ImGui::Dummy(ImVec2(0.0f, 15.0f));

        ImGui::Text("Impurity settings");
        ImGui::SliderScalar("Region", ImGuiDataType_Double, &sp.region_size, &region_bounds.x, &region_bounds.y, "%.2e");
        ImGui::SliderScalar("Extends", ImGuiDataType_Double, &sp.region_extends, &extends_bounds.x, &extends_bounds.y, "%.2e");
        ImGui::SliderScalar("Density", ImGuiDataType_Double, &sp.impurity_density, &density_bounds.x, &density_bounds.y, "%.2e");
        ImGui::SliderScalar("Radius", ImGuiDataType_Double, &sp.impurity_radius, &radius_bounds.x, &radius_bounds.y, "%.2e");

        if (ImGui::Button("Random seed")) {
            sp.impurity_seed += 1;
            compute_requested = true;
        }

        ImGui::Dummy(ImVec2(0.0f, 20.0f));

        compute_requested |= ImGui::Button("Compute"); ImGui::SameLine();
        ImGui::Checkbox("Sync immediately", &sync_immediate);
    }
    ImGui::End();

    return parameters_changed || compute_requested;
}


void TexturesView::ShowView(const std::vector<uint32_t> & textures)
{
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("Output");

    auto textureId = textures[0];
    auto viewportSize = ImGui::GetContentRegionAvail();
    float min_dim = viewportSize.x < viewportSize.y ? viewportSize.x : viewportSize.y;
    ImGui::Image((void*)textureId, ImVec2(min_dim, min_dim));

    ImGui::PopStyleVar();
    ImGui::End();
}
