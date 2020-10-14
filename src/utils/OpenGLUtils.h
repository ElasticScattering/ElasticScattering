#pragma once

#include <string>
#include <GL/glew.h>

class OpenGLUtils {
	std::string ReadShaderFile(const char* shader_file);

public:
	void Init(GLuint& vbo, GLuint& vao, GLuint& shader_program);
};

