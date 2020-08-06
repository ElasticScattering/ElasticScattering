#define DOCTEST_CONFIG_IMPLEMENT

#include <windows.h>

#include "ElasticScattering.h"
#include "utils/Test.h"
#include "utils/ErrorMacros.h"
#include "doctest.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>


enum class Mode {
    LIFETIME,
    STATS,
    AVG_DISTANCE,
    CONDUCTIVITY
};

typedef struct
{
    bool run_tests;
    int num_iterations;
    Mode mode;
    bool show_info;
} InitParameters;

void ProcessInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

std::string ReadShaderFile(const char *shader_file) 
{
    std::ifstream file(shader_file);
    std::stringstream sstream;
    sstream << file.rdbuf();

    std::string contents = sstream.str();
    return contents;
}

void ParseArgs(int argc, char** argv, InitParameters* p_init) {
    p_init->num_iterations = 1;
    p_init->mode = Mode::LIFETIME;
    p_init->show_info = true;
    p_init->run_tests = false;

    if (argc == 2) {
        
        std::string w;
        w.assign(argv[1], strlen(argv[1]));

        p_init->run_tests = w == "test";
        if (p_init->run_tests)
            std::cout << "Testing! " << std::endl;
    }
    return;
    
    if (argc != 4) {
        std::cout << "Usage: ElasticScattering [num iterations] [lifetime | distance | stats | conductivity] [show | no-show]" << std::endl;
        std::cout << "Usage: ElasticScattering [test]" << std::endl;
        exit(0);
    }

    p_init->num_iterations = atoi(argv[1]);

    std::unordered_map<std::string, Mode> modes
    {
        {"lifetime",     Mode::LIFETIME},
        {"distance",     Mode::AVG_DISTANCE},
        {"stats",        Mode::STATS},
        {"conductivity", Mode::CONDUCTIVITY}
    };
    std::string key;
    key.assign(argv[2], strlen(argv[2]));
    auto iterator = modes.find(key);
    ERR_FAIL_COND_MSG(iterator == modes.end(), "Couldn't understand second command line argument.");
    p_init->mode = iterator->second;

    p_init->show_info = strcmp(argv[3], "show");
}

int main(int argc, char **argv)
{
    InitParameters init;
    ParseArgs(argc, argv, &init);

    if (init.run_tests)
    {
        doctest::Context context;
        context.applyCommandLine(argc, argv);
        context.setOption("no-breaks", true);
        int res = context.run();

        return 0;
    }

    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    int width = 900;
    int height = 900;
    window = glfwCreateWindow(width, height, "Elastic Scattering", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    RECT rect;
    GetClientRect(GetDesktopWindow(), &rect);
    glfwSetWindowPos(window, rect.right / 2 - width / 2, rect.bottom / 2 - height / 2);

    glfwMakeContextCurrent(window);
    glViewport(0, 0, width, height);

    SimulationParameters sp;
    sp.region_size      = 1e-6;
    sp.particle_count   = 10'000; //100'000'000;
    sp.particle_row_count = sqrt(sp.particle_count);
    sp.particle_speed   = 7e5;
    sp.particle_mass    = 5 * M0;
    sp.impurity_count   = 100;
    sp.impurity_radius  = 1.5e-8;
    sp.impurity_radius_sq = sp.impurity_radius * sp.impurity_radius;
    sp.alpha            = PI / 4.0;
    sp.phi = 0;// sp.alpha - 1e-10;
    sp.magnetic_field   = 20;
    sp.angular_speed    = E * sp.magnetic_field / sp.particle_mass;
    sp.tau              = 1e-12;
    
    std::cout << "\n\n+---------------------------------------------------+" << std::endl;
    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "Start region size: (" << sp.region_size << ", " << sp.region_size << ")" << std::endl;
    std::cout << "Particles:         " << sp.particle_count << std::endl;
    std::cout << "Particle speed:    " << sp.particle_speed << std::endl;
    std::cout << "Particle mass:     " << sp.particle_mass << std::endl;
    std::cout << "Impurities:        " << sp.impurity_count << std::endl;
    std::cout << "Impurity radius:   " << sp.impurity_radius << std::endl;
    std::cout << "Alpha:             " << sp.alpha << std::endl;
    std::cout << "Phi:               " << sp.phi << std::endl;
    std::cout << "Magnetic field:    " << sp.magnetic_field << std::endl;
    std::cout << "Angular speed:     " << sp.angular_speed << std::endl;
    std::cout << "Tau:               " << sp.tau << std::endl;
    std::cout << "-----------------------------------------------------" << std::endl;

    ERR_FAIL_COND_MSG(pow(sp.particle_row_count, 2) != sp.particle_count, "Particles couldn't be placed in a square grid");
    ERR_FAIL_COND_MSG(sp.alpha > (PI / 4.0), "Alpha should not be greater than pi/4.");
    ERR_FAIL_COND_MSG(sp.alpha <= 0, "Alpha should be positive.");
    ERR_FAIL_COND_MSG(sp.angular_speed < 0, "Angular speed (w) should be positive");
    ERR_FAIL_COND_MSG(sp.magnetic_field < 0, "Magnetic field strength (B) should be positive");

    ElasticScattering* es = new CPUElasticScattering();
    es->Init(sp);
    es->Compute();
    auto pixels = es->GetPixels();

    std::cout << "\n\n+---------------------------------------------------+" << std::endl;

    GLenum error = glewInit();
    if (error != GLEW_OK) return EXIT_FAILURE;

    // Shaders
    GLint success;

    std::string source = ReadShaderFile("shader.vs");
    char *vsource = new char[source.length() + 1];
    std::copy_n(source.c_str(), source.length() + 1, vsource);

    source = ReadShaderFile("shader.fs");
    char* fsource = new char[source.length() + 1];
    std::copy_n(source.c_str(), source.length() + 1, fsource);

    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, &vsource, nullptr);
    glCompileShader(vert_shader);
    
    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[2048];
        glGetShaderInfoLog(vert_shader, 2048, nullptr, infoLog);
        std::cout << "Failed to compile vertex shader. Info:\n" << infoLog << std::endl;
    }

    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(frag_shader, 1, &fsource, nullptr);
    glCompileShader(frag_shader);

    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[2048];
        glGetShaderInfoLog(frag_shader, 2048, nullptr, infoLog);
        std::cout << "Failed to compile fragment shader. Info:\n" << infoLog << std::endl;
    }

    GLuint shader_program = glCreateProgram();
    glAttachShader(shader_program, vert_shader);
    glAttachShader(shader_program, frag_shader);
    glLinkProgram(shader_program);
    glUseProgram(shader_program);
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader);

    // Vertex data 
    float vertices[] =
    {
        // Position,        Tex coord
        -0.9f, -0.9f, 0.0f, 0.0f, 0.0f,
         0.9f, -0.9f, 0.0f, 1.0f, 0.0f,
         0.9f,  0.9f, 0.0f, 1.0f, 1.0f,
        -0.9f,  0.9f, 0.0f, 0.0f, 1.0f
    };

    GLuint vbo, vao;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    
    glBindVertexArray(vao);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    
    auto stride = 5 * sizeof(float);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Texture
    int dim = sqrt(pixels.size()/3);
    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dim, dim, 0, GL_RGB, GL_FLOAT, pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
#if 1
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
#else
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
#endif
    
    glUseProgram(shader_program);
    glUniform1i(glGetUniformLocation(shader_program, "texture1"), 0);

    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, tex);

        glUseProgram(shader_program);
        glBindVertexArray(vao);
        glDrawArrays(GL_QUADS, 0, 4);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    glDeleteProgram(shader_program);
    glfwTerminate();
    return 0;
}