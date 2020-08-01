#include <windows.h>

#include "ElasticScattering.h"

#include<string>
#include<iostream>
#include<sstream>
#include<fstream>

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glfw3.h>

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


int main(void)
{
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

    ElasticScattering* es = new ElasticScattering();
    es->Init(0, nullptr);

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

    double* data = es->GetData();
    int dim = 100; // sqrt(sizeof(data) / sizeof(*data));

    int length = dim * dim;
    double itau = 1.0 / 1e-12; // es->result_max_time;
    float *pixels = new float[length * 3];
    int j = 0;
    for (int i = 0; i < length; i++)
    {
        double k = data[i] * itau;
        if (k == 0) {
            pixels[j] = 1.0f;
            pixels[j + 1] = 0.0f;
            pixels[j + 2] = 0.0f;
        }
        else {
            pixels[j] = k;
            pixels[j + 1] = k;
            pixels[j + 2] = k;
        }
        j += 3;
    }

    GLuint tex;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dim, dim, 0, GL_RGB, GL_FLOAT, pixels);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

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

        /*
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, dim, dim, 0, GL_RGBA, GL_DOUBLE, data);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 2, 2, 0, GL_RGB, GL_FLOAT, pixels);

        glBindTexture(GL_TEXTURE_2D, texID);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

       glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f, -1.0f);
        glTexCoord2f(1.0f, 1.0f); glVertex2f(1.0f, -1.0f);
        glTexCoord2f(1.0f, 0.0f); glVertex2f(1.0f, 1.0f);
        glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, 1.0f);
        glEnd();
        */

        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    glDeleteProgram(shader_program);

    es->Cleanup();

    glfwTerminate();
    return 0;
}