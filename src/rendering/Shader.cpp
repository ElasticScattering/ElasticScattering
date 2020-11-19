#pragma once

#include "Shader.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <GL/glew.h>
#include <GL/glfw3.h>

void Shader::SetInt(const std::string& name, int value)
{
    GLint location = glGetUniformLocation(id, name.c_str());
    glUniform1i(location, value);
}

char* Shader::ReadFile(const std::string& shader_file)
{
    std::ifstream file(shader_file);
    std::stringstream sstream;
    sstream << file.rdbuf();

    std::string contents = sstream.str();

    char* vsource = new char[contents.length() + 1];
    std::copy_n(contents.c_str(), contents.length() + 1, vsource);

    return vsource;
}

Shader::Shader(const std::string& vertexSrc, const std::string& fragmentSrc)
{
    auto vsource = ReadFile("src/shaders/shader.vs");
    auto fsource = ReadFile("src/shaders/shader.fs");

    GLuint vert_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vert_shader, 1, &vsource, nullptr);
    glCompileShader(vert_shader);

    GLint success;

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

    id = glCreateProgram();
    glAttachShader(id, vert_shader);
    glAttachShader(id, frag_shader);
    
    // Could get info log here.
    glLinkProgram(id);

    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader);
}

Shader::~Shader()
{
    glDeleteProgram(id);
}

void Shader::Bind() const
{
    glUseProgram(id);
}

void Shader::Unbind() const
{
    glUseProgram(0);
}
