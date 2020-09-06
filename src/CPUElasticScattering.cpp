#include "ElasticScattering.h"
#include "OpenGLUtils.h"

double CPUElasticScattering::Compute(const SimulationParameters &p_sp)
{
    if (!PrepareCompute(p_sp)) return 0;

    int limit = sp.dim - 1;

    for (int j = 0; j < limit; j++) {
        for (int i = 0; i < limit; i++)
        {
            v2 pos(i, j);
            pos = pos * (sp.region_size / (double)(sp.dim - 2));

            double particle_result = (sp.mode == MODE_DIR_LIFETIME) ? single_lifetime(pos, sp.phi, &sp, impurities) : phi_lifetime(pos, &sp, impurities);
            main_buffer[j * sp.dim + i] = particle_result * GetWeight2D(i, j, sp.dim);
        }
    }

    double result = ComputeResult(main_buffer);
    MakeTexture();

    return result;
}

bool CPUElasticScattering::PrepareCompute(const SimulationParameters &p_sp) {
    if (!first_run && !AnythingChanged(p_sp)) return false;

    bool impurities_changed = ImpuritySettingsChanged(p_sp);
    bool work_size_changed  = (sp.dim != p_sp.dim);

    sp = p_sp;
    CompleteSimulationParameters();
        
    if (first_run || impurities_changed)
        GenerateImpurities(sp);

    if (first_run || work_size_changed) {
        main_buffer.clear();
        main_buffer.resize(particle_count, 0);
    }

    first_run = false;

    return true;
}

void CPUElasticScattering::MakeTexture() 
{
    pixels.clear();
    pixels.resize(particle_count * 4L);

    int j = 0;
    for (int i = 0; i < particle_count; i++)
    {
        float k = GetColor(main_buffer[i], sp.tau, sp.mode);

        if (sp.mode != MODE_SIGMA_XY) {
            pixels[j]      = k;
            pixels[j + 1L] = k;
            pixels[j + 2L] = k;
        }
        else {
            if (k < 0.0) {
                pixels[j] = 0;
                pixels[j + 1L] = 0;
                pixels[j + 2L] = -k;
            } else {
                pixels[j] = k;
                pixels[j + 1L] = 0;
                pixels[j + 2L] = 0;
            }
        }
        pixels[j + 3L] = 1.0;
        j += 4;
    }

    glGenTextures(1, &ogl.tex);
    glBindTexture(GL_TEXTURE_2D, ogl.tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, sp.dim, sp.dim, 0, GL_RGBA, GL_FLOAT, pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glUseProgram(ogl.shader_program);
    glUniform1i(glGetUniformLocation(ogl.shader_program, "texture1"), 0);
}

CPUElasticScattering::CPUElasticScattering()
{
    OpenGLUtils o;
    o.Init(ogl.vbo, ogl.vao, ogl.shader_program);
}